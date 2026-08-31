#include "region.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace pld2 {

namespace {

constexpr int64_t kMaxPos = INT64_MAX - 1;

bool file_exists(const std::string& path) {
    std::ifstream f(path.c_str());
    return f.good();
}

// Parses one "chr:start-end" / "chr:start" / "chr" token into `r`.
// start/end are 1-based inclusive. Returns false + err on malformed input.
bool parse_region_token(const std::string& tok, Region& r, std::string& err) {
    const size_t colon = tok.find(':');
    if (colon == std::string::npos) {
        r.chrom = tok;
        if (r.chrom.empty()) {
            err = "-L: empty region";
            return false;
        }
        r.start = 1;
        r.end = kMaxPos;
        r.label = r.chrom;
        return true;
    }
    r.chrom = tok.substr(0, colon);
    if (r.chrom.empty()) {
        err = "-L: empty chromosome in \"" + tok + "\"";
        return false;
    }
    const std::string coord = tok.substr(colon + 1);
    const size_t dash = coord.find('-');
    auto parse_num = [&](const std::string& s, int64_t& v) -> bool {
        if (s.empty()) return false;
        for (char c : s)
            if (!std::isdigit(static_cast<unsigned char>(c))) return false;
        v = strtoll(s.c_str(), nullptr, 10);
        return v >= 1;
    };
    int64_t a, b;
    if (dash == std::string::npos) {
        if (!parse_num(coord, a)) {
            err = "-L: bad coordinate in \"" + tok + "\"";
            return false;
        }
        b = a;  // single position -> 1 bp region
    } else {
        if (!parse_num(coord.substr(0, dash), a) ||
            !parse_num(coord.substr(dash + 1), b)) {
            err = "-L: bad coordinate in \"" + tok + "\"";
            return false;
        }
    }
    if (b < a) {
        err = "-L: start > end in \"" + tok + "\"";
        return false;
    }
    r.start = a;
    r.end = b;
    r.label = r.chrom + ":" + std::to_string(a) + "-" + std::to_string(b);
    return true;
}

bool parse_bed(const std::string& path, std::vector<Region>& out,
               std::string& err) {
    std::ifstream f(path.c_str());
    if (!f) {
        err = "-L: cannot open BED file \"" + path + "\"";
        return false;
    }
    std::string line;
    int line_no = 0;
    while (std::getline(f, line)) {
        line_no++;
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        // chrom<TAB>start<TAB>end  (0-based half-open); extra columns ignored.
        size_t p0 = line.find('\t');
        if (p0 == std::string::npos) p0 = line.find(' ');
        if (p0 == std::string::npos) {
            std::cerr << "warning: -L BED line " << line_no
                      << " has no coordinates, skipped\n";
            continue;
        }
        std::string chrom = line.substr(0, p0);
        size_t p1 = line.find('\t', p0 + 1);
        if (p1 == std::string::npos) p1 = line.find(' ', p0 + 1);
        if (p1 == std::string::npos) {
            std::cerr << "warning: -L BED line " << line_no
                      << " has no end column, skipped\n";
            continue;
        }
        auto num = [](const std::string& s) -> int64_t {
            if (s.empty()) return -1;
            for (char c : s)
                if (!std::isdigit(static_cast<unsigned char>(c))) return -1;
            return strtoll(s.c_str(), nullptr, 10);
        };
        const int64_t start0 = num(line.substr(p0 + 1, p1 - p0 - 1));
        size_t p2 = line.find_first_of("\t ", p1 + 1);
        const std::string end_str =
            (p2 == std::string::npos) ? line.substr(p1 + 1)
                                      : line.substr(p1 + 1, p2 - p1 - 1);
        const int64_t end0 = num(end_str);
        if (start0 < 0 || end0 < 0 || start0 >= end0) {
            std::cerr << "warning: -L BED line " << line_no
                      << " invalid interval, skipped\n";
            continue;
        }
        Region r;
        r.chrom = chrom;
        r.start = start0 + 1;  // 0-based half-open -> 1-based inclusive
        r.end = end0;
        r.label = r.chrom + ":" + std::to_string(r.start) + "-" +
                  std::to_string(r.end);
        out.push_back(std::move(r));
    }
    if (out.empty()) {
        err = "-L: BED file \"" + path + "\" contains no usable regions";
        return false;
    }
    return true;
}

}  // namespace

bool load_regions(const std::string& spec, std::vector<Region>& out,
                  std::string& err) {
    out.clear();
    if (spec.empty()) return true;

    std::vector<Region> parsed;
    if (file_exists(spec)) {
        if (!parse_bed(spec, parsed, err)) return false;
    } else {
        size_t start = 0;
        while (start <= spec.size()) {
            size_t comma = spec.find(',', start);
            if (comma == std::string::npos) comma = spec.size();
            std::string tok = spec.substr(start, comma - start);
            Region r;
            if (!parse_region_token(tok, r, err)) return false;
            parsed.push_back(std::move(r));
            if (comma == spec.size()) break;
            start = comma + 1;
        }
    }

    // Canonical deterministic order used by both the fallback and fast paths.
    std::sort(parsed.begin(), parsed.end(),
              [](const Region& x, const Region& y) { return x.label < y.label; });

    // Drop exact duplicates (same chrom + coordinates -> same label, adjacent
    // after the sort). Without this, a repeated region in the -L argument or a
    // duplicated BED row would be an independent accumulation unit on the fast
    // path (one task per region) while the fallback folds both into a single
    // label-keyed buffer -- double-counting every pair vs the fallback.
    auto last = std::unique(parsed.begin(), parsed.end(),
                            [](const Region& x, const Region& y) {
                                return x.label == y.label;
                            });
    parsed.erase(last, parsed.end());

    out = std::move(parsed);
    return true;
}

void RegionIndex::add(int64_t start, int64_t end, const std::string& label) {
    starts_.push_back(start);
    ends_.push_back(end);
    labels_.push_back(&label);
}

const std::string* RegionIndex::find(int64_t pos) const {
    // Last region with start <= pos (upper_bound - 1).
    auto it = std::upper_bound(starts_.begin(), starts_.end(), pos);
    if (it == starts_.begin()) return nullptr;
    const size_t idx = static_cast<size_t>(it - starts_.begin()) - 1;
    if (pos > ends_[idx]) return nullptr;
    return labels_[idx];
}

RegionSet::RegionSet(const std::vector<Region>& regions) {    // Group by chrom; keep insertion order within a chrom (then binary-search
    // on starts; overlaps resolve to the largest start).
    std::vector<const Region*> order;
    order.reserve(regions.size());
    for (const Region& r : regions) order.push_back(&r);
    std::sort(order.begin(), order.end(),
              [](const Region* x, const Region* y) {
                  if (x->chrom != y->chrom) return x->chrom < y->chrom;
                  return x->start < y->start;
              });
    std::string cur;
    for (const Region* r : order) {
        if (r->chrom != cur) {
            cur = r->chrom;
            by_chrom_.emplace(cur, RegionIndex{});
        }
        by_chrom_.find(cur)->second.add(r->start, r->end, r->label);
    }
}

std::string region_query_string(const Region& r) {
    if (r.start == 1 && r.end == kMaxPos) return r.chrom;
    return r.chrom + ":" + std::to_string(r.start) + "-" + std::to_string(r.end);
}

}  // namespace pld2