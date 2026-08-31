#include "output.h"

#include <cstdio>
#include <cstring>
#include <iostream>

#include <htslib/bgzf.h>

namespace pld2 {

int write_stat(const std::string& path, const std::vector<Bin>& bins,
               int max_dist_bp, int out_type) {
    BGZF* fp = bgzf_open(path.c_str(), "w");
    if (!fp) {
        std::cerr << "open OUT File error: " << path << std::endl;
        return 0;
    }

    const char* header = "#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\n";
    if (bgzf_write(fp, header, static_cast<int>(strlen(header))) < 0) {
        std::cerr << "write OUT File error: " << path << std::endl;
        bgzf_close(fp);
        std::remove(path.c_str());  // don't leave a truncated file behind
        return 0;
    }

    // Original OUTStatFile: real D' columns only when 2<=TF<=5 (in practice
    // only OutType 2 reaches it with real values; 4/5 use the histogram path).
    const bool real_d = (out_type >= 2 && out_type <= 5);
    char buf[256];
    for (int ii = 1; ii <= max_dist_bp; ++ii) {
        const Bin& b = bins[static_cast<size_t>(ii)];
        if (b.count == 0) continue;
        const double meanRR = b.sumRR / static_cast<double>(b.count);
        int len;
        if (real_d) {
            const double meanD = b.sumD / static_cast<double>(b.count);
            len = snprintf(buf, sizeof(buf), "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%lld\n",
                           ii, meanRR, meanD, b.sumRR, b.sumD,
                           static_cast<long long>(b.count));
        } else {
            len = snprintf(buf, sizeof(buf), "%d\t%.4f\tNA\t%.4f\tNA\t%lld\n",
                           ii, meanRR, b.sumRR, static_cast<long long>(b.count));
        }
        if (bgzf_write(fp, buf, len) < 0) {
            std::cerr << "write OUT File error: " << path << std::endl;
            bgzf_close(fp);
            std::remove(path.c_str());  // don't leave a truncated file behind
            return 0;
        }
    }

    bgzf_close(fp);
    return 1;
}

int write_hist_stat(const std::string& path, int max_dist_bp,
                    const std::vector<int>& hist_rr,
                    const std::vector<int>& hist_d, size_t nhaplo,
                    int out_type) {
    BGZF* fp = bgzf_open(path.c_str(), "w");
    if (!fp) {
        std::cerr << "open OUT File error: " << path << std::endl;
        return 0;
    }

    const char* header =
        (out_type == 5) ? "#Dist\tR^2\tR^2_count\n"
                        : "#Dist\tR^2\tR^2_count\tD'\tD_count\n";
    char buf[512];
    if (bgzf_write(fp, header, static_cast<int>(strlen(header))) < 0) {
        bgzf_close(fp);
        std::remove(path.c_str());  // don't leave a truncated file behind
        return 0;
    }
    int len = snprintf(buf, sizeof(buf), "#2SampleSize\t%zu\n", nhaplo);
    if (bgzf_write(fp, buf, len) < 0) {
        std::cerr << "write OUT File error: " << path << std::endl;
        bgzf_close(fp);
        std::remove(path.c_str());  // don't leave a truncated file behind
        return 0;
    }

    // Rows: for dist 1..max_dist_bp, bucket j=0..100 with a non-zero count.
    // The D' value column reuses the r^2 bucket value (j/100) -- the original
    // prints RR_double there; kept for byte-compatibility.
    const bool has_d = (out_type == 4);

    // Sparse scan: only dists with at least one non-zero bucket are visited,
    // in ascending dist order (so the output rows are unchanged).
    std::vector<int> touched;
    touched.reserve(256);
    for (int i = 0; i < max_dist_bp; ++i) {
        const size_t base = static_cast<size_t>(i) * 101;
        for (int j = 0; j < 101; ++j) {
            if (hist_rr[base + j] != 0 || (has_d && hist_d[base + j] != 0)) {
                touched.push_back(i);
                break;
            }
        }
    }

    for (int i : touched) {
        for (int j = 0; j < 101; ++j) {
            const size_t idx = static_cast<size_t>(i) * 101 + j;
            const int rr_cnt = hist_rr[idx];
            const int d_cnt = has_d ? hist_d[idx] : 0;
            if (rr_cnt == 0 && d_cnt == 0) continue;
            const double rr_double = static_cast<double>(j) / 100.0;
            if (has_d) {
                len = snprintf(buf, sizeof(buf), "%d\t%.2f\t%d\t%.2f\t%d\n",
                               i + 1, rr_double, rr_cnt, rr_double, d_cnt);
            } else {
                len = snprintf(buf, sizeof(buf), "%d\t%.2f\t%d\n", i + 1,
                               rr_double, rr_cnt);
            }
            if (bgzf_write(fp, buf, len) < 0) {
                std::cerr << "write OUT File error: " << path << std::endl;
                bgzf_close(fp);
                std::remove(path.c_str());  // don't leave a truncated file behind
                return 0;
            }
        }
    }

    bgzf_close(fp);
    return 1;
}

}  // namespace pld2
