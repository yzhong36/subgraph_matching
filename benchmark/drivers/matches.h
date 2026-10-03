// Optional match output of the batch drivers, for verification runs (off unless an output path is given).
// One line per distinct node set among a query's matches: the query index, then the set's target node indices, sorted,
// separated by spaces.  Node indices are positions in the target file, counted from 0.
#pragma once
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>

inline void write_node_set(FILE* f, int query, const std::vector<int>& nodes) {
    fprintf(f, "%d", query);
    for (int x : nodes) fprintf(f, " %d", x);
    fputc('\n', f);
}

struct NodeSets {
    struct Hash {
        size_t operator()(const std::vector<int>& v) const {
            size_t h = 1469598103934665603ULL;
            for (int x : v) h ^= (size_t)x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_set<std::vector<int>, Hash> sets;
    template <class It> void add(It first, It last) {
        std::vector<int> v;
        for (; first != last; ++first) v.push_back((int)*first);
        std::sort(v.begin(), v.end());
        sets.insert(std::move(v));
    }
    void write(FILE* f, int query) const { for (const auto& s : sets) write_node_set(f, query, s); }
};

// worker w's part of the matches file, and the parent's merge of the parts
inline std::string matches_part(const std::string& path, int w) { return path + "." + std::to_string(w); }
inline void merge_matches(const std::string& path, int workers) {
    std::ofstream out(path);
    for (int w = 0; w < workers; w++) {
        std::ifstream in(matches_part(path, w));
        if (in.peek() != std::ifstream::traits_type::eof()) out << in.rdbuf();
        std::remove(matches_part(path, w).c_str());
    }
}
