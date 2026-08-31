#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "subpop.h"
#include "types.h"

namespace pld2 {

// One genotype site: physical position + encoded haplotypes (0 = best, 1 =
// second, 2 = missing/unknown), in input order.
struct GenoSite {
    int64_t pos = 0;
    std::vector<uint8_t> vals;
};

// Unpacked per-contig accumulation buffer from the genotype reader (sites in
// file order; packing/sorting happens in merge_geno_buffers).
struct GenoBuf {
    std::string name;
    std::vector<GenoSite> sites;
};

// Reads a PopLDdecay `-InGenotype` file into unpacked per-contig buffers,
// replicating the original two paths:
//   - non-SubPop: FilterGenotype.h FilterGeno (Cut3base always on in the main
//     tool) + GetBestBase() encoding. With -OutFilterSNP the passing lines are
//     also streamed verbatim to `{prefix}.genotype.filter.gz` (plain gzip,
//     byte-compatible with the original ogzstream output).
//   - SubPop: Read_SubPopGenotype_IN() inline filter (miss/het denominators
//     over the subset, MAF = SeD/(SeD+Max), biallelic-only) + encode (het ->
//     best/sed pseudo-haplotypes, missing -> 2). No filter.gz here.
// The two paths intentionally differ (MAF definition, missing encoding), as in
// the original. Prints the original stdout/stderr messages. With SubPop specs
// each group is resolved against the genotype #CHROM header and packed
// independently (one output vector per group); with no specs a single
// all-samples vector is produced. Returns 1 on success, 0 on error (header
// missing, unopenable files, subset < 3, ...).
int read_genotype(const std::string& path, const Options& opt,
                  const std::string& prefix,
                  const std::vector<SubGroupSpec>& specs,
                  std::vector<std::vector<GenoBuf>>& out);

// Merges the unpacked genotype buffers into `out` (typically packed ContigData
// from the VCF reader; may be empty). Matching contig names are combined (sites
// merged, stable-sorted by pos, dedupe keep-first with the VCF site winning a
// tie); new names are appended. Per-chr n_haplotypes = the first (smallest-pos)
// site's haplotype count, replicating the original "Asize = first site" LD
// loop; other sites are padded with best (0) / truncated.
void merge_geno_buffers(std::vector<GenoBuf>& gbufs, std::vector<ContigData>& out);

}  // namespace pld2