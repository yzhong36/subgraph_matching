// ri_batch: RI 3.6 (Bonnici et al., BMC Bioinformatics 2013), the authors' code, monomorphism, for many queries on one
// target.  Target and every query read once, before the clock; W forked workers take one query at a time from a shared
// counter.  Per query, inside the clock, what ri3.cpp's match() does for "mono gfu": the matching machine
// (MaMaConstrFirst) and match() with the authors' EmptyMatchListener, which counts the matches.
//   ri_batch target.gfu queries.txt workers out.tsv [limit_s [matches_out]]
//   matches_out (verification runs, not for timing): the matches go to a listener that keeps their distinct node sets
//   (in place of ri3.cpp's ConsoleMatchListener, which prints them), written as in matches.h
//   out.tsv: query index, seconds, matches (mappings), timed out (0/1); stdout: wall seconds.
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
#include "rows.h"

using namespace rilib;

// verification runs only: counts as EmptyMatchListener does and keeps each match's node set
class SetListener : public MatchListener {
public:
    NodeSets sets;
    virtual void match(int n, int* qIDs, int* rIDs) { matchcount++; sets.add(rIDs, rIDs + n); }
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
    // every worker's output files are opened before the clock starts, so file-system latency is not timed
    std::vector<FILE*> fos(workers), fms(workers, nullptr);
    for (int w = 0; w < workers; w++) {
        fos[w] = fopen((out + "." + std::to_string(w)).c_str(), "w");
        if (!matches.empty()) fms[w] = fopen(matches_part(matches, w).c_str(), "w");
    }
    auto t0 = Clock::now();
    for (int w = 0; w < workers; w++) {
        if (fork() == 0) {
            FILE* fo = fos[w];
            FILE* fm = fms[w];
            signal(SIGALRM, on_alarm);
            StringAttrComparator nodeCmp; DefaultAttrComparator edgeCmp;
            std::string rows;
            for (;;) {
                int j = next->fetch_add(1);
                if (j >= nq) break;
                long mappings = 0; int timed_out = 0;
                auto s0 = std::chrono::steady_clock::now();
                SetListener* S = fm ? new SetListener() : nullptr;
                MatchListener* L = S ? (MatchListener*)S : new EmptyMatchListener();
                if (sigsetjmp(jb, 1) == 0) {
                    alarm(limit);
                    MaMaConstrFirst* mama = new MaMaConstrFirst(*queries[j]);
                    mama->build(*queries[j]);
                    long st, tc, mc;
                    match(*target, *queries[j], *mama, *L, MT_MONO, nodeCmp, edgeCmp, &st, &tc, &mc);
                    alarm(0);
                    mappings = L->matchcount;
                    if (S) S->sets.write(fm, j);
                    delete mama; delete L;
                } else {
                    timed_out = 1; mappings = L->matchcount;                            // L and the solver leak
                }
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
                { char b[192]; snprintf(b, sizeof b, "%d\t%.6f\t%ld\t%d\n", j, sec, mappings, timed_out); rows += b; }
            }
            finish_worker(fo, rows, t0);
            if (fm) fclose(fm);
            _exit(0);
        }
    }
    while (wait(nullptr) > 0) {}
    double wall = std::chrono::duration<double>(Clock::now() - t0).count();
    for (FILE* f : fos) fclose(f);
    for (FILE* f : fms) if (f) fclose(f);
    wall = merge_rows(out, workers, wall);
    if (!matches.empty()) merge_matches(matches, workers);
    std::cout << wall << "\n";
    return 0;
}
