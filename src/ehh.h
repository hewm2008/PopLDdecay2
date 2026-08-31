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

// Parses the -EHH spec by splitting on ':'. The original bails out with the
// "Para [-EHH] should be [chr:Site]" message from the CLI when the colon is
// missing; the caller is expected to have validated that already, so here we
// only split (empty tokens skipped, like the original split()).
bool parse_ehh(const std::string& spec, EhhParams& out);

// Runs the EHH decay for the region [start_site-max_dist_bp, start_site+
// max_dist_bp] on the matching contig, writing `{prefix}.ehh.gz` (plain gzip,
// byte-compatible with the original ogzstream output), then the R plot script
// / Rscript run and the `{prefix}.stat.gz` removal. `tf` mirrors the original
// -OutType value controlling whether the temp R script is deleted.
//
// Returns 1 on the full success path and 0 on the bail-out paths (missing
// chromosome, <8 no-missing SNPs) -- mirroring EHH_Region_LDDecay. On the
// bail-out paths an empty .ehh.gz is still created and the .stat.gz file is
// left in place, exactly like the original.
int run_ehh(const std::vector<ContigData>& contigs, const EhhParams& p,
            int max_dist_bp, int tf, const std::string& prefix);

}  // namespace pld2
