// vf3_batch: VF3 / VF3L (vf3lib, MiviaLab; Carletti et al., TPAMI 2018), the authors' code, for many queries on one
// target.  Monomorphism (-e), undirected (-u).  Target and every query read once, before the clock; W forked workers take
// one query at a time from a shared counter.  Per query, inside the clock, what main.cpp does for one run (-r 0):
// FastCheck, node classification of both graphs, VF3 node sorting, and FindAllMatchings.
// Guard (not in main.cpp): a pattern using a label the target lacks has no match and is answered 0 without the search,
// whose node classifier would give that label a class beyond the target's and index past its arrays (labels.h; the
// lookups are inside the query's time).
//   build: benchmark/setup_tools.sh (-DVF3 or -DVF3L)
//   vf3_batch target.vft queries.txt workers out.tsv [limit_s [matches_out]]
//   matches_out (verification runs, not for timing): also write every query's distinct matched node sets (matches.h)
//   out.tsv: query index, seconds, solutions (mappings), timed out (0/1); stdout: wall seconds.
#include <atomic>
#include <chrono>
#include <csetjmp>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include "VFLib.h"
#include "labels.h"
#include "matches.h"
#include "rows.h"
#include "Options.hpp"
using namespace vflib;

static sigjmp_buf jb;
static void on_alarm(int) { siglongjmp(jb, 1); }

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: vf3_batch target.vft queries.txt workers out.tsv [limit_s [matches_out]]\n"; return 1; }
    int workers = atoi(argv[3]);
    std::string out = argv[4];
    int limit_s = argc > 5 ? atoi(argv[5]) : 300;
    std::string matches = argc > 6 ? argv[6] : "";
    Options opt; opt.undirected = true; opt.edgeInduced = true; opt.format = "vf";
    opt.storeSolutions = !matches.empty();
    std::ifstream tin(argv[1]);
    ARGLoader<data_t, Empty>* tl = CreateLoader<data_t, Empty>(opt, tin);
    ARGraph<data_t, Empty> targ(tl);
    std::vector<ARGraph<data_t, Empty>*> pats;
    std::vector<std::string> qpaths;
    { std::ifstream f(argv[2]); std::string s;
      while (std::getline(f, s)) if (!s.empty()) { qpaths.push_back(s); std::ifstream pin(s); pats.push_back(new ARGraph<data_t, Empty>(CreateLoader<data_t, Empty>(opt, pin))); } }
    int nq = pats.size();
    LabelGuard guard(argv[1], qpaths, "vft");
    std::atomic<int>* next = (std::atomic<int>*)mmap(nullptr, sizeof(std::atomic<int>), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    new (next) std::atomic<int>(0);
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
            MatchingEngine<state_t>* me = CreateMatchingEngine(opt);
            std::string rows;
            for (;;) {
                int j = next->fetch_add(1);
                if (j >= nq) break;
                size_t sols = 0; int timed_out = 0;
                auto s0 = std::chrono::steady_clock::now();
                if (sigsetjmp(jb, 1) == 0) {
                    alarm(limit_s);
                    ARGraph<data_t, Empty>& patt = *pats[j];
                    me->ResetSolutionCounter();
                    FastCheck<data_t, data_t, Empty, Empty> check(&patt, &targ);
                    if (guard.possible(j) && check.CheckSubgraphIsomorphism()) {
                        NodeClassifier<data_t, Empty> classifier(&targ);
                        NodeClassifier<data_t, Empty> classifier2(&patt, classifier);
                        std::vector<uint32_t> class_patt = classifier2.GetClasses();
                        std::vector<uint32_t> class_targ = classifier.GetClasses();
                        uint32_t classes_count = classifier.CountClasses();
                        VF3NodeSorter<data_t, Empty, SubIsoNodeProbability<data_t, Empty>> sorter(&targ);
                        std::vector<nodeID_t> sorted = sorter.SortNodes(&patt);
                        state_t st(&patt, &targ, class_patt.data(), class_targ.data(), classes_count, sorted.data(), opt.edgeInduced);
                        me->FindAllMatchings(st);
                        sols = me->GetSolutionsCount();
                        if (fm) {
                            std::vector<MatchingSolution> found; me->GetSolutions(found);
                            NodeSets ns;
                            for (auto& sol : found) { std::vector<int> t; for (auto& pr : sol) t.push_back((int)pr.second); ns.add(t.begin(), t.end()); }
                            ns.write(fm, j);
                        }
                    }
                    alarm(0);
                } else {
                    timed_out = 1; sols = me->GetSolutionsCount();                      // the state leaks
                }
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
                { char b[192]; snprintf(b, sizeof b, "%d\t%.6f\t%zu\t%d\n", j, sec, sols, timed_out); rows += b; }
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
