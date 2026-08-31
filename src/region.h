#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace pld2 {

// One -L region. Coordinates are 1-based INCLUSIVE on both ends (VCF/htslib
// convention). A whole-chromosome region has start=1, end=INT64_MAX-1.
struct Region {
    std::string chrom;  // contig name (as it appears in the VCF)
    int64_t start = 1;  // 1-based inclusive
    int64_t end = 0;    // 1-based inclusive
    std::string label;  // canonical display key: "chr" or "chr:start-end"
};

// Parses the -L argument into `out`. The argument is either:
//   - a comma-separated list of regions:  "chr", "chr:start", "chr:start-end"
//     (start/end 1-based inclusive), or
//   - a path to an existing BED file (0-based half-open [start,end) per line;
//     extra columns ignored; '#'/blank lines skipped).
// Detection: if the string names an existing file it is read as BED, otherwise
// parsed as a region list. Output is sorted by `label` so the fallback and fast
// paths process regions in the same deterministic order. Returns false (with an
// error message) on malformed input.
bool load_regions(const std::string& spec, std::vector<Region>& out,
                  std::string& err);

// Per-chrom region membership lookup (sorted starts + binary search).
class RegionIndex {
  public:
    // Build from the regions belonging to this chrom.
    RegionIndex() = default;
    void add(int64_t start, int64_t end, const std::string& label);

    // 1-based inclusive position -> containing region label, or nullptr.
    const std::string* find(int64_t pos) const;

  private:
    std::vector<int64_t> starts_;  // sorted ascending (1-based starts)
    std::vector<int64_t> ends_;    // parallel: inclusive ends
    std::vector<const std::string*> labels_;
};

// Fast membership lookup over all regions. For overlapping regions on the same
// chrom the region with the largest start wins (deterministic).
class RegionSet {
  public:
    explicit RegionSet(const std::vector<Region>& regions);

    // 1-based inclusive position on `chrom` -> containing region label or "".
    const std::string& find(const std::string& chrom, int64_t pos) const {
        auto it = by_chrom_.find(chrom);
        if (it == by_chrom_.end()) return empty_;
        const std::string* l = it->second.find(pos);
        return l ? *l : empty_;
    }

    // Length-bounded lookup for a contig name inside a record buffer (no
    // temporary std::string): the map uses a transparent comparator so a
    // std::string_view compares directly against the std::string keys.
    const std::string& find(const char* chrom, size_t len, int64_t pos) const {
        auto it = by_chrom_.find(std::string_view(chrom, len));
        if (it == by_chrom_.end()) return empty_;
        const std::string* l = it->second.find(pos);
        return l ? *l : empty_;
    }

  private:
    // Transparent comparator (std::less<>) so find() accepts string_view keys.
    std::map<std::string, RegionIndex, std::less<>> by_chrom_;
    std::string empty_;
};

// Tabix query string for a region: "chr" for a whole contig, otherwise the
// 1-based inclusive "chr:start-end" (matches htslib's region-string convention).
std::string region_query_string(const Region& r);

}  // namespace pld2