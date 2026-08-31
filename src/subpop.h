#pragma once

#include <set>
#include <string>
#include <vector>

namespace pld2 {

// Resolved -SubPop sample subset (one per population).
struct SubPop {
    bool active = false;           // true when the -SubPop group resolved ok
    std::string label;             // output suffix ("" = legacy single group)
    std::vector<int> sample_cols;  // VCF token indices (>=9) of the subset
                                   // samples, in VCF header order
    int n_samples = 0;             // == sample_cols.size() (NumberSubGroup)
};

// A -SubPop group before it is resolved against a header: a label and the
// requested sample names.
struct SubGroupSpec {
    std::string label;   // "" for the legacy single-column single-file group
    std::set<std::string> names;
};

// Parses the -SubPop file list into groups. Each file is either
// single-column (one group per file; label = stripped basename) or
// two-column `sample<TAB>groupid` (one group per distinct groupid; label =
// that id). When exactly one single-column file is given the label is left
// empty (legacy mode: output name unchanged). '#'/blank lines and trailing
// '\r' are skipped, like the original sample-list reader. Returns false if a
// file cannot be opened or the column widths are inconsistent.
bool parse_subpop_lists(const std::vector<std::string>& paths,
                        std::vector<SubGroupSpec>& out);

// Resolves group specs against the VCF header (read once here), replicating
// the original Read_SubPopVCF_IN* prelude: whitespace-separated names per
// line matched to the #CHROM columns, subset = the VCF header order of
// matched samples. Prints the original per-group messages ("the Number of
// subPop samples[found in VCF] is N", not-found/repeat warnings). Returns
// false on error: header missing/broken, or any group <3 samples.
bool resolve_subpops(const std::string& vcf_path,
                     const std::vector<SubGroupSpec>& specs,
                     std::vector<SubPop>& out);

// Resolves group specs against an already-parsed #CHROM header token list
// (genotype input), printing the per-group messages with the given found-msg
// ("[found in VCF]" / "[found in Genotype]") and header-kind in warnings
// ("in the VCF Header" / "in the Genotype Header"). `sample_start` is the
// first sample column (9 for VCF, 2 for genotype). Returns false if any group
// <3 samples.
bool resolve_subpop_header(const std::vector<std::string>& header,
                           const std::vector<SubGroupSpec>& specs,
                           const char* found_msg, const char* header_kind,
                           int sample_start, std::vector<SubPop>& out);

}  // namespace pld2