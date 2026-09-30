// ri_batch: RI 3.6 monomorphism (the classical non-induced subgraph matching) for many queries on one target.
// The target and every query are read once, before the clock; then W forked workers take one query at a time
// from a shared counter (the accounting of seeded exact's multiprocessing pool, chunksize 1).
// Per query: the RI matching order is built and the search run inside the clock; every match's reference node
// set is sorted and kept in a hash set, so the count is of distinct node sets, like the ground truth.
//   ri_batch target.gfu queries.txt workers out.tsv [limit_s [matches_out]]
//   matches_out (verification runs, not for timing): also write every query's distinct matched node sets (matches.h)
//   out.tsv: query index, seconds, mappings, distinct node sets, timed out (0/1); stdout: wall seconds.
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <atomic>
#include <csetjmp>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

#include "c_textdb_driver.h"
#include "timer.h"
#include "AttributeComparator.h"
#include "AttributeDeconstructor.h"
#include "Graph.h"
#include "MatchingMachine.h"
#include "MaMaConstrFirst.h"
#include "Match.h"
#include "Solver.h"
#include "SubGISolver.h"
#include "matches.h"

using namespace rilib;

struct VecHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t h = 1469598103934665603ULL;
        for (int x : v) { h ^= (size_t)x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); }
        return h;
    }
};

class SetListener : public MatchListener {
public:
    std::unordered_set<std::vector<int>, VecHash> sets;
    virtual void match(int n, int* qIDs, int* rIDs) {
        matchcount++;
        std::vector<int> v(rIDs, rIDs + n);
        std::sort(v.begin(), v.end());
        sets.insert(std::move(v));
    }
};

static sigjmp_buf jb;
static void on_alarm(int) { siglongjmp(jb, 1); }

int main(int argc, char* argv[]) {
    if (argc < 5) { std::cerr << "usage: ri_batch target.gfu queries.txt workers out.tsv [limit_s [matches_out]]\n"; return 1; }
    int workers = atoi(argv[3]);
    int limit = argc > 5 ? atoi(argv[5]) : 300;
    std::string matches = argc > 6 ? argv[6] : "";
    Graph* target = new Graph();
    if (read_graph(argv[1], target, GFT_GFU) != 0) { std::cerr << "cannot read target\n"; return 1; }
    std::vector<std::string> qpaths;
    { std::ifstream f(argv[2]); std::string s; while (std::getline(f, s)) if (!s.empty()) qpaths.push_back(s); }
    std::vector<Graph*> queries;
    for (auto& p : qpaths) { Graph* q = new Graph(); if (read_graph(p.c_str(), q, GFT_GFU) != 0) { std::cerr << "bad " << p << "\n"; return 1; } queries.push_back(q); }
    int nq = queries.size();
    std::atomic<int>* next = (std::atomic<int>*)mmap(nullptr, sizeof(std::atomic<int>), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    new (next) std::atomic<int>(0);
    std::string out = argv[4];
    auto t0 = std::chrono::steady_clock::now();
    for (int w = 0; w < workers; w++) {
        if (fork() == 0) {
            FILE* fo = fopen((out + "." + std::to_string(w)).c_str(), "w");
            FILE* fm = matches.empty() ? nullptr : fopen(matches_part(matches, w).c_str(), "w");
            signal(SIGALRM, on_alarm);
            StringAttrComparator nodeCmp; DefaultAttrComparator edgeCmp;
            for (;;) {
                int j = next->fetch_add(1);
                if (j >= nq) break;
                long mappings = 0, sets = 0; int timed_out = 0;
                auto s0 = std::chrono::steady_clock::now();
                SetListener* L = new SetListener();
                if (sigsetjmp(jb, 1) == 0) {
                    alarm(limit);
                    MaMaConstrFirst* mama = new MaMaConstrFirst(*queries[j]);
                    mama->build(*queries[j]);
                    long st, tc, mc;
                    match(*target, *queries[j], *mama, *L, MT_MONO, nodeCmp, edgeCmp, &st, &tc, &mc);
                    alarm(0);
                    mappings = L->matchcount; sets = L->sets.size();
                    if (fm) for (const auto& s : L->sets) write_node_set(fm, j, s);
                    delete mama; delete L;
                } else {
                    timed_out = 1; mappings = L->matchcount; sets = L->sets.size();     // L and the solver leak
                }
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
                fprintf(fo, "%d\t%.6f\t%ld\t%ld\t%d\n", j, sec, mappings, sets, timed_out);
                fflush(fo);
            }
            fclose(fo);
            if (fm) fclose(fm);
            _exit(0);
        }
    }
    while (wait(nullptr) > 0) {}
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream fo(out);
    for (int w = 0; w < workers; w++) {
        std::ifstream f(out + "." + std::to_string(w)); std::string s;
        while (std::getline(f, s)) fo << s << "\n";
        std::remove((out + "." + std::to_string(w)).c_str());
    }
    if (!matches.empty()) merge_matches(matches, workers);
    std::cout << wall << "\n";
    return 0;
}
