// Label pre-check shared by the batch drivers.  A pattern that uses a label the target lacks has no match, so the query is
// answered 0 without calling the tool.  DAF and VF3 need this (DAF reads past its label table for such a label, VF3's node
// classifier would index past the target's classes); RI and FASTiso get it too, so every tool is treated alike.  The label
// lists are read with the query files, before the clock; only the set lookups (microseconds) run inside a query's time.
#pragma once
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// The node labels of a graph file.  fmt: "gfu" ('#g', n, then n label lines), "vft" (n, then n lines "id label"),
// "daf" ('v id label' lines).
inline std::vector<long> file_labels(const std::string& path, const std::string& fmt) {
    std::vector<long> out;
    std::ifstream f(path);
    std::string line;
    if (fmt == "daf") {
        while (std::getline(f, line))
            if (line.size() > 1 && line[0] == 'v') { std::istringstream ss(line.substr(1)); long id, l; ss >> id >> l; out.push_back(l); }
        return out;
    }
    if (fmt == "gfu") std::getline(f, line);                            // "#g"
    std::getline(f, line);
    long n = std::stol(line);
    for (long i = 0; i < n && std::getline(f, line); i++) {
        std::istringstream ss(line);
        long a, b;
        ss >> a;
        if (fmt == "vft") { ss >> b; out.push_back(b); } else out.push_back(a);
    }
    return out;
}

struct LabelGuard {
    std::set<long> target;
    std::vector<std::vector<long>> queries;
    LabelGuard(const std::string& target_path, const std::vector<std::string>& query_paths, const std::string& fmt) {
        for (long l : file_labels(target_path, fmt)) target.insert(l);
        for (const auto& q : query_paths) queries.push_back(file_labels(q, fmt));
    }
    bool possible(int j) const {
        for (long l : queries[j]) if (!target.count(l)) return false;
        return true;
    }
};
