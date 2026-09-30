// Result rows of the batch drivers.  A worker keeps its rows in memory and writes them once after its last query, with the
// moment that query finished ("#done <seconds since the clock started>").  The reported wall time is when the last query
// finished, so it contains neither a file-system write per query nor process teardown.
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>

using Clock = std::chrono::steady_clock;      // system-wide monotonic: a worker's reading is comparable with the parent's

inline void finish_worker(FILE* fo, const std::string& rows, Clock::time_point t0) {
    double done = std::chrono::duration<double>(Clock::now() - t0).count();
    fputs(rows.c_str(), fo);
    fprintf(fo, "#done\t%.6f\n", done);
    fclose(fo);
}

// Merge the workers' row files into `out`; returns when the last query finished, or `fallback` if a worker did not report.
inline double merge_rows(const std::string& out, int workers, double fallback) {
    std::ofstream fo(out);
    double last = 0;
    int reported = 0;
    for (int w = 0; w < workers; w++) {
        std::string part = out + "." + std::to_string(w), s;
        std::ifstream f(part);
        while (std::getline(f, s)) {
            if (s.rfind("#done", 0) == 0) { last = std::max(last, std::stod(s.substr(6))); reported++; }
            else fo << s << "\n";
        }
        std::remove(part.c_str());
    }
    return reported == workers ? last : fallback;
}
