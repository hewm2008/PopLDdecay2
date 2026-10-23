#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "types.h"

namespace pld2 {

// `-EHH chr:site` core-site spec (replicates the original `chr:Site`).
struct EhhParams {
    std::string chr;
    int64_t start_site = 0;
};

// Parses the -EHH spec by splitting on ':' (empty tokens skipped, like the
// original split()) and requires the site to be a whole number >= 0.
bool parse_ehh(const std::string& spec, EhhParams& out);

// Runs the EHH decay for the region [start_site-max_dist_bp, start_site+
// max_dist_bp] on the matching contig, writing `{prefix}.ehh.gz` (plain gzip,
// byte-compatible with the original ogzstream output), then the R plot script
// / Rscript run and the `{prefix}.stat.gz` removal. `tf` mirrors the original
// -OutType value controlling whether the temp R script is deleted.
//
// Returns 1 on the full success path and 0 on the bail-out paths (missing
// chromosome, <8 no-missing SNPs, or a write failure) -- mirroring
// EHH_Region_LDDecay. On the two original bail-outs an empty .ehh.gz is still
// created and .stat.gz is left in place, exactly like the original; on a
// write failure the .ehh.gz is removed.
int run_ehh(const std::vector<ContigData>& contigs, const EhhParams& p,
            int max_dist_bp, int tf, const std::string& prefix);

}  // namespace pld2
