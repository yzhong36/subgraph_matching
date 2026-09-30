// daf_batch: DAF (Han, Kim, Gu, Park & Han, SIGMOD 2019), the authors' code, for many queries on one data graph.
// The data graph is read and processed once and every query file read, before the clock; W forked workers take one
// query at a time from a shared counter.  Per query, inside the clock, what main/main.cc times: query processing
// against the data graph, the DAG, the candidate space and the backtracking (all matches).
// Two guards, not in main.cc, both for queries with no embedding:
//  - a query label above the data graph's largest label is answered 0 without DAF.  DAF's label table
//    (DataGraph::GetTransferredLabel) marks the labels up to the largest that the data graph lacks, but is read past its
//    end for a larger one: the authors' program then segfaults or reports matches that do not exist.  The labels are
//    read with the query files, before the clock; the comparison is inside the query's time.
//  - when QueryGraph::LoadAndProcessGraph reports a label the data graph lacks, or the candidate space cannot be built,
//    the count is 0 and the partly built objects are left alone, since destroying them crashes (free(): invalid pointer).
//   daf_batch data.daf queries.txt workers out.tsv [limit_s]
//   (no match output: DAF counts the matches of degree-1 query nodes combinatorially instead of enumerating them)
//   out.tsv: query index, seconds, matches, timed out (0/1); stdout: wall seconds.
#include <atomic>
#include <chrono>
#include <csetjmp>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <sys/mman.h>
#include <sys/wait.h>
#include <sstream>
#include <unistd.h>
#include "rows.h"
#include "include/backtrack.h"
#include "include/candidate_space.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

static sigjmp_buf jb;
static void on_alarm(int) { siglongjmp(jb, 1); }

// the largest vertex label of a DAF graph file ("v id label" lines)
static long max_label(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    long top = -1;
    while (std::getline(f, line))
        if (line.size() > 1 && line[0] == 'v') { std::istringstream ss(line.substr(1)); long id, l; ss >> id >> l; top = std::max(top, l); }
    return top;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: daf_batch data.daf queries.txt workers out.tsv [limit_s]\n"; return 1; }
    int workers = atoi(argv[3]);
    std::string out = argv[4];
    int limit_s = argc > 5 ? atoi(argv[5]) : 300;
    daf::DataGraph data(argv[1]);
    data.LoadAndProcessGraph();
    std::vector<daf::QueryGraph*> queries;
    std::vector<std::string> qpaths;
    { std::ifstream f(argv[2]); std::string s;
      while (std::getline(f, s)) if (!s.empty()) { queries.push_back(new daf::QueryGraph(s)); qpaths.push_back(s); } }
    int nq = queries.size();
    long data_top = max_label(argv[1]);
    std::vector<long> query_top;
    for (auto& q : qpaths) query_top.push_back(max_label(q));
    std::atomic<int>* next = (std::atomic<int>*)mmap(nullptr, sizeof(std::atomic<int>), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    new (next) std::atomic<int>(0);
    // every worker's output file is opened before the clock starts, so file-system latency is not timed
    std::vector<FILE*> fos(workers);
    for (int w = 0; w < workers; w++) fos[w] = fopen((out + "." + std::to_string(w)).c_str(), "w");
    auto t0 = Clock::now();
    for (int w = 0; w < workers; w++) {
        if (fork() == 0) {
            FILE* fo = fos[w];
            if (!freopen("/dev/null", "w", stdout)) return 1;
            signal(SIGALRM, on_alarm);
            std::string rows;
            for (;;) {
                int j = next->fetch_add(1);
                if (j >= nq) break;
                uint64_t count = 0; int timed_out = 0;
                auto s0 = std::chrono::steady_clock::now();
                if (sigsetjmp(jb, 1) == 0) {
                    alarm(limit_s);
                    daf::QueryGraph& query = *queries[j];
                    if (query_top[j] <= data_top && query.LoadAndProcessGraph(data)) {
                        auto* dag = new daf::DAG(data, query);
                        dag->BuildDAG();
                        auto* cs = new daf::CandidateSpace(data, query, *dag);
                        if (cs->BuildCS()) {
                            auto* bt = new daf::Backtrack(data, query, *cs);
                            count = bt->FindMatches(std::numeric_limits<uint64_t>::max());
                            delete bt; delete cs; delete dag;
                        }                                                   // else: no embedding; objects left alone
                    }
                    alarm(0);
                } else timed_out = 1;
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
                { char b[192]; snprintf(b, sizeof b, "%d\t%.6f\t%lu\t%d\n", j, sec, (unsigned long)count, timed_out); rows += b; }
            }
            finish_worker(fo, rows, t0);
            _exit(0);
        }
    }
    while (wait(nullptr) > 0) {}
    double wall = std::chrono::duration<double>(Clock::now() - t0).count();
    for (FILE* f : fos) fclose(f);
    wall = merge_rows(out, workers, wall);
    std::cout << wall << "\n";
    return 0;
}
