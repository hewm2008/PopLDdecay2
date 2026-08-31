#pragma once

#include <string>
#include <vector>

#include "subpop.h"
#include "types.h"

namespace pld2 {

struct AccumSpec;

// Reads a VCF (plain / gzip / bgzip) via htslib and packs filtered sites into
// per-contig SoA arrays. Replicates the original Read_VCF_IN / Read_VCF_IN_Phase
// filter + encode semantics exactly (see doc/PITFALLS.md). With SubPop groups,
// filters + packing run per group over that group's subset samples (like the
// original Read_SubPopVCF_IN*); with no groups a single all-samples set is
// produced. `outs` gets one ContigData vector per group. Returns 1 on success,
// 0 on error.
int read_vcf(const std::string& path, const Options& opt,
             std::vector<std::vector<ContigData>>& outs, ReadStats& stats,
             const std::vector<SubPop>& subs);

// Fast path (M3a): reads + filters + packs + computes each chromosome in a
// dedicated worker thread, using the tabix index (.tbi/.csi) for per-chromosome
// random access. Read of one chromosome overlaps the LD compute of another.
// Requires: index present, opt.threads > 1, >=1 chromosome in the header's
// ##contig lines (single-chromosome inputs are admitted too and position-shard
// when -T > 1). Tasks are built from the index's own seqnames (tbx_seqnames),
// so contigs absent from ##contig lines are still processed; names that would
// break a tabix region query (empty or containing ':') or a failed tbx_itr_querys
// cause a decline (returns false) so no data is silently dropped. Returns false
// (caller falls back to read_vcf+compute_bins) when those conditions are not
// met. On success fills `outs` (one AccumSpec per group) and `stats`; the
// caller must set outs[].out_type and size outs to the number of groups first.
bool compute_bins_indexed(const std::string& path, const Options& opt,
                          int max_dist_bp, std::vector<AccumSpec>& outs,
                          ReadStats& stats, const std::vector<SubPop>& subs);

}  // namespace pld2