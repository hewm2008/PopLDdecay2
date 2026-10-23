#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "ldwriter.h"
#include "types.h"

namespace pld2 {

// Per-contig / per-run LD accumulation context. Depending on out_type the
// per-pair values feed distance bins (1,2,3,6,7,8), value histograms
// (4,5), and/or per-pair .LD.gz rows (3,6,7,8); OutType 0 computes no pair
// at all (its .stat.gz is a header only). `nhaplo` = 2*N (sub-population
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
// sites ascending by pos, pairs (i<j) ascending, and the per-pair math
// op-for-op the original cal_RR_* formula (see ldkernel.cpp) -> byte-identical.
//
// `core_n` (default: all sites) limits which sites act as the LEFT partner of a
// pair; the right partners still run over the whole buffer. That is how the
// fastpath shards positions: each shard owns the pairs whose left site falls
// in its core window (with max_dist_bp of context on the right), so every pair
// is counted exactly once. `part_label` (default: cd.name) is the LdWriter
// part label; a fastpath shard passes a sortable composite label (chr + shard
// index) so finish() concatenates the parts in globally ascending (chr, p1)
// order.
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
// enforced by RunTest/06.Verify/golden_diff.sh -T tests). `out.nhaplo` must be
// set beforehand.
void compute_bins(const std::vector<ContigData>& contigs, int max_dist_bp,
                  int nthreads, AccumSpec& out);

// Prepares one AccumSpec for a run: the distance bins, the OutType 4/5
// histograms and, for the pairwise modes, the .LD.gz writer.
void init_accum(AccumSpec& spec, int out_type, size_t max_dist_bp);

// Caps the number of accumulation workers in the OutType 4/5 histogram mode:
// each worker holds a (MaxDist+1)*101 int histogram (~121 MB at MaxDist=300 kb).
// The sums are integers, so they merge in any order and the output stays
// byte-identical.
void cap_histogram_workers(size_t& nworkers, int out_type);

// Merges one worker's partial result into the final AccumSpec. The integer
// counts are exact; the double sums are added in worker order (<=1 ULP from
// the sequential order, which the %.4f output does not show).
void merge_partial(AccumSpec& out, const AccumSpec& partial);

// Progress reporting, on stderr (see the stdout/stderr rule in main.cpp).
// The counters are shared by the worker threads, so the numbers are the totals
// over the whole run, not per thread; the contig name shown is the one the
// reporting thread is working on. A line is printed at most every 5 minutes.
void ld_progress_reset();

}  // namespace pld2
