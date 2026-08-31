#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "ldwriter.h"
#include "types.h"

namespace pld2 {

// Per-contig / per-run LD accumulation context. Depending on out_type the
// per-pair values feed distance bins (0,1,2,3,6,7,8), value histograms
// (4,5), and/or per-pair .LD.gz rows (3,6,7,8). `nhaplo` = 2*N (sub-population
// adjusted; used for the OutType 4/5 "#2SampleSize" header line).
struct AccumSpec {
    int out_type = 1;
    size_t nhaplo = 0;         // 2*N_eff
    std::vector<Bin> bins;     // distance bins, indexed by bp (1..max_dist_bp)
    std::vector<int> hist_rr;  // OutType 4/5: (max_dist_bp+1)*101 histogram
    std::vector<int> hist_d;   // OutType 4:   (max_dist_bp+1)*101 histogram
    std::unique_ptr<LdWriter> ld;  // OutType 3/6/7/8
};

// Accumulates LD for one contig into `spec` (bins and/or histograms and/or
// .LD.gz rows), replicating the original PairWiseComV1 per-contig loop exactly:
// sites ascending by pos, pairs (i<j) ascending. The per-pair math dispatches
// on spec.out_type and replicates the corresponding original cal_RR_* formula
// op-for-op (see ldkernel.cpp), so output is byte-identical to the original.
//
// `core_n` (default: all sites) limits which sites act as the LEFT partner of a
// pair (sites [0, core_n)); right partners still run over the whole buffer. This
// supports fastpath position sharding, where each shard owns pairs whose left
// site falls in its core window (overlapping max_dist_bp of context on the
// right), so every pair is counted exactly once. `part_label` (default: cd.name)
// is the LdWriter part label used for .LD.gz ordering; fastpath shards pass a
// sortable composite label (chr + shard index) so finish() concatenates parts in
// globally ascending (chr, p1) order.
void accumulate_contig(const ContigData& cd, int max_dist_bp, AccumSpec& spec,
                       size_t core_n = SIZE_MAX,
                       const std::string& part_label = {});

// Computes LD for OutType 0-8 into `out`. With nthreads<=1 (or a single
// non-empty contig) accumulation is strictly sequential over contigs in input
// order (byte-identical to the original). With nthreads>1, non-empty contigs
// are distributed to worker threads (greedy largest-first to the least-loaded
// worker); each worker accumulates its own partial bins/histograms/.LD parts
// and the results are merged. Integer counts exact; double sums ordered per
// worker (may differ from sequential by <=1 ULP; %.4f output byte-identical,
// enforced by golden_diff.sh -T tests). `out.nhaplo` must be set beforehand.
void compute_bins(const std::vector<ContigData>& contigs, int max_dist_bp,
                  int nthreads, AccumSpec& out);

}  // namespace pld2
