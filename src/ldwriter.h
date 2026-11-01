#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <zlib.h>

namespace pld2 {

// Per-pair .LD.gz output (OutType 3/6/7/8), replicating the original tool's
// `{prefix}.LD.gz` (plain gzip via ogzstream). Rows are per (chr, site1, site2)
// in ascending order. Each instance owns a private temp directory and buffers
// each contig in its own gzip-member temp file, so the parallel fast path can
// accumulate rows concurrently. finish() writes the header as its own member
// and then concatenates the temp members in lexicographic label order (the
// original's std::map chr order), then cleans up.
class LdWriter {
  public:
    explicit LdWriter(int out_type);
    ~LdWriter();

    LdWriter(const LdWriter&) = delete;
    LdWriter& operator=(const LdWriter&) = delete;

    // Switch row output to the temp part for `label` (closes the previous).
    // Call once per contig before emitting its rows.
    void begin_contig(const std::string& label);

    // OutType 3.
    void emit_rr(const std::string& chr, int64_t p1, int64_t p2, double rr,
                 int64_t dis);
    // OutType 6.
    void emit_d_rr(const std::string& chr, int64_t p1, int64_t p2, double d,
                   double rr, int64_t dis);
    // OutType 7.
    void emit_d_lod_rr(const std::string& chr, int64_t p1, int64_t p2, double d,
                       double lod, double rr, int64_t dis);
    // OutType 8 (CI columns printed fixed %.2f of low_i/high_i over 100).
    void emit_d_lod_rr_ci(const std::string& chr, int64_t p1, int64_t p2,
                          double d, double lod, double rr, int low_i,
                          int high_i, int64_t dis);

    // Moves `other`'s buffered parts into this writer (fast-path merge).
    void absorb(LdWriter& other);

    // Streams all buffered parts into `{prefix}.LD.gz` (plain gzip made of
    // concatenated gzip members, one per part, in label order). Returns true
    // on success; on a lost part it writes nothing and returns false.
    bool finish(const std::string& prefix);

  private:
    void close_part();

    int out_type_;
    std::string tmpdir_;
    gzFile cur_ = nullptr;
    int part_seq_ = 0;
    // Set when a part could not be opened, written or closed, or when a part
    // was lost while merging. finish() then refuses to publish a table that
    // would silently be missing rows.
    bool failed_ = false;
    struct Part {
        std::string label;
        std::string path;
    };
    std::vector<Part> parts_;
};

}  // namespace pld2
