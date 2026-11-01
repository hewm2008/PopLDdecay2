#pragma once

#include <string>
#include <vector>

#include "types.h"

namespace pld2 {

// Writes `{path}` (a .stat.gz path) for the distance-bin OutTypes
// (0,1,2,3,6,7,8), replicating the original OUTStatFile byte-for-byte: the
// header, a row per bin 1..max_dist_bp with count>0, and real Mean_D'/Sum_D'
// columns only where the original has them (the others print "NA").
// Returns 1 on success, 0 on error.
int write_stat(const std::string& path, const std::vector<Bin>& bins,
               int max_dist_bp, int out_type);

// Writes the OutType 4/5 histogram .stat.gz (replicating the original
// PairWiseComNewOUT_A): "#Dist\tR^2\tR^2_count[D'\tD_count]" + a
// "#2SampleSize\t{2*N}" line, then rows for dist 1..max_dist_bp and value
// bucket j=0..100 (printed as j/100 with %.2f, D' value column reusing the
// r^2 bucket value -- an original bug kept for byte-compat). `hist` holds
// (max_dist_bp+1)*101 ints. Returns 1 on success, 0 on error.
int write_hist_stat(const std::string& path, int max_dist_bp,
                    const std::vector<int>& hist_rr,
                    const std::vector<int>& hist_d, size_t nhaplo,
                    int out_type);

}  // namespace pld2
