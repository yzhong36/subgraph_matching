// fastiso_batch: FASTiso (Coti, Reinharz & Agbeto, bioRxiv 2025), the authors' code, for many queries on one target.
// Monomorphism, RI's gfu format read as undirected (as '-A monomorphism -F gfu -G -S count-solution').  Target and every
// query parsed once, before the clock; W forked workers take one query at a time from a shared counter.  Per query,
// inside the clock, what fast_iso.cpp's main does after parsing: initial domains (Filter<CheckMonomorphism>), node
// ordering (NodeSorter, TRI_ORDER 0) and SearchMatching counting all solutions.
//   fastiso_batch target.gfu queries.txt workers out.tsv [limit_s [matches_out]]
//   matches_out (verification runs, not for timing): also write every query's distinct matched node sets (matches.h);
//   FASTiso then runs with '-S print-solution' and keeps its solutions in memory instead of counting
//   out.tsv: query index, seconds, solutions (mappings), timed out (0/1); stdout: wall seconds.
#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <unordered_map> 
#include <ParseArgs.hpp>
#include <Parser/Parser.hpp>
#include <Parser/ParseFile.hpp>
#include <Parser/ParseCSVFile.hpp>
#include <Parser/ParseLADFile.hpp>
#include <Parser/ParseGRFFile.hpp>
#include <Parser/ParseGFUFile.hpp>
#include <Graph.hpp>
#include <Checker/CheckIso.hpp>
#include <Checker/CheckIsomorphism.hpp>
#include <Checker/CheckSubIsomorphism.hpp>
#include <Checker/CheckMonomorphism.hpp>
#include <SearchMatching.hpp>
#include <Filter.hpp>
#include <utils.hpp>
#include <sbitset.h>
#include <NodeSorter/NodeSorted.hpp>
#include <NodeSorter/RINodeSorted.hpp>
#include <NodeSorter/VF3NodeSorted.hpp>
#include <State/State.hpp>
#include <State/NodeInfo.hpp>
#include <time_utils.hpp>
#include <unistd.h>

using namespace std;

namespace fastiso::internal
{
    #ifdef WIN32
    LARGE_INTEGER _freq;
    #else
    struct timezone _tz;
    #endif
}

using namespace fastiso::internal;


#include <atomic>
#include <csetjmp>
#include <csignal>
#include <fstream>
#include <sys/mman.h>
#include <sys/wait.h>
#include "matches.h"

static sigjmp_buf jb;
static void on_alarm(int) { siglongjmp(jb, 1); }

static size_t run_query(Graph& pattern, Graph& target, const CompareEmptyEdgeAttr& cmp,
                        std::vector<std::vector<State::Monomorphism::node_id_nt>>* keep = nullptr) {
    sbitset* data = new sbitset[pattern.size()];
    vector<node_id_t> data_size(pattern.size());
    Filter<Checker::CheckMonomorphism> filter;
    bool ok = filter.filterNode(data, data_size, pattern, target, ParseArgs::getFILTER_NODE());
    if (!ok) { delete[] data; return 0; }
    vector<size_t> dm_result(pattern.size(), 0);
    vector<int> node_parent(pattern.size(), -1);
    NodeSorter::NodeSorter nodeSort;
    nodeSort.sortNodes(pattern, node_parent, data_size, dm_result, false);
    std::vector<node_id_t>& node_order = data_size;
    SearchMatching<State::Monomorphism::SubState, State::Monomorphism::NodeInfo, 4> sm(ParseArgs::getAll_solution(), ParseArgs::getSOLUTION(), cmp);
    sm.solutions = keep;
    sm.searchMatching(data, node_order, node_parent, pattern, target, dm_result);
    size_t count = sm.num_result;
    delete[] data;
    return count;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: fastiso_batch target.gfu queries.txt workers out.tsv [limit_s [matches_out]]\n"; return 1; }
    int workers = atoi(argv[3]);
    std::string out = argv[4];
    int limit_s = argc > 5 ? atoi(argv[5]) : 300;
    std::string matches = argc > 6 ? argv[6] : "";
    std::vector<std::string> qpaths;
    { std::ifstream f(argv[2]); std::string s; while (std::getline(f, s)) if (!s.empty()) qpaths.push_back(s); }
    // the command line's options, so every static setting the solver reads is what the CLI would use
    std::string a0 = "FASTiso", a1 = "-A", a2 = "monomorphism", a3 = "-F", a4 = "gfu", a5 = "-G", a6 = "-S", a7 = matches.empty() ? "count-solution" : "print-solution";
    std::string a8 = qpaths[0], a9 = argv[1];
    char* fake[] = {a0.data(), a1.data(), a2.data(), a3.data(), a4.data(), a5.data(), a6.data(), a7.data(), a8.data(), a9.data(), nullptr};
    ParseArgs::fast_args_parse(10, fake);
    Graph target(0);
    { auto* pf = new Parser::ParseGFUFile(target, true, false); Parser::Parser::parse(*pf, argv[1]); delete pf; }
    std::vector<Graph*> pats;
    for (auto& q : qpaths) { Graph* g = new Graph(0); auto* pf = new Parser::ParseGFUFile(*g, true, false); Parser::Parser::parse(*pf, q); delete pf; pats.push_back(g); }
    const CompareEmptyEdgeAttr* cmp = new const CompareEmptyEdgeAttr();
    int nq = pats.size();
    std::atomic<int>* next = (std::atomic<int>*)mmap(nullptr, sizeof(std::atomic<int>), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    new (next) std::atomic<int>(0);
    auto t0 = std::chrono::steady_clock::now();
    for (int w = 0; w < workers; w++) {
        if (fork() == 0) {
            FILE* fo = fopen((out + "." + std::to_string(w)).c_str(), "w");
            FILE* fm = matches.empty() ? nullptr : fopen(matches_part(matches, w).c_str(), "w");
            if (!freopen("/dev/null", "w", stdout)) return 1;
            signal(SIGALRM, on_alarm);
            for (;;) {
                int j = next->fetch_add(1);
                if (j >= nq) break;
                size_t count = 0; int timed_out = 0;
                auto s0 = std::chrono::steady_clock::now();
                std::vector<std::vector<State::Monomorphism::node_id_nt>> found;
                if (sigsetjmp(jb, 1) == 0) {
                    alarm(limit_s); count = run_query(*pats[j], target, *cmp, fm ? &found : nullptr); alarm(0);
                    if (fm) { NodeSets ns; for (auto& m : found) ns.add(m.begin(), m.end()); ns.write(fm, j); }
                } else timed_out = 1;
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
                fprintf(fo, "%d\t%.6f\t%zu\t%d\n", j, sec, count, timed_out);
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
