#include "ldwriter.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

#include <zlib.h>

#include "rowfmt.h"
#include "textutil.h"

namespace pld2 {

namespace {

// Single-quotes a path for use in a system() shell command (TMPDIR comes from
// the environment and may contain spaces / shell metacharacters).

const char* header_for(int out_type) {
    switch (out_type) {
        case 3: return "#chr\tSite1\tSite2\tr^2\tDist\n";
        case 6: return "#chr\tSite1\tSite2\tD'\tr^2\tDist\n";
        case 7: return "#chr\tSite1\tSite2\tD'\tLOD\tr^2\tDist\n";
        case 8: return "#chr\tSite1\tSite2\tD'\tLOD\tr^2\tCIlow\tCIhi\tDist\n";
        default: return "";
    }
}

// Removes a directory and everything below it (temp cleanup).
void remove_tree(const std::string& path) {
    // Only ever removes our own mkdtemp dirs (no symlink traversal).
    std::string cmd = "rm -rf -- " + shell_quote(path);
    if (std::system(cmd.c_str()) != 0) {
        std::cerr << "warning: could not clean temp dir " << path << std::endl;
    }
}

// Compresses `data` into a standalone gzip member (deflateInit2 windowBits
// 15+16 emits the full gzip format: magic header, deflate stream, CRC+ISIZE
// trailer). Returns false on a zlib error.
bool gzip_member(const std::string& data, std::string& out) {
    z_stream zs{};
    if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK) {
        return false;
    }
    const uLong bound = deflateBound(&zs, static_cast<uLong>(data.size()));
    std::string buf(bound, '\0');
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data.data()));
    zs.avail_in = static_cast<uInt>(data.size());
    zs.next_out = reinterpret_cast<Bytef*>(&buf[0]);
    zs.avail_out = static_cast<uInt>(bound);
    const int rc = deflate(&zs, Z_FINISH);
    deflateEnd(&zs);
    if (rc != Z_STREAM_END) return false;
    buf.resize(zs.total_out);
    out = std::move(buf);
    return true;
}

}  // namespace

LdWriter::LdWriter(int out_type) : out_type_(out_type) {
    std::string base = "/tmp";
    if (const char* td = std::getenv("TMPDIR")) {
        if (*td) base = td;
    }
    std::string tmpl = base + "/pld2.ld.XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    char* made = mkdtemp(buf.data());
    if (made) tmpdir_ = made;
}

LdWriter::~LdWriter() {
    close_part();
    if (!tmpdir_.empty()) remove_tree(tmpdir_);
}

void LdWriter::begin_contig(const std::string& label) {
    close_part();
    if (tmpdir_.empty()) {
        std::cerr << "error: could not create LD temp directory (TMPDIR="
                  << (std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp")
                  << "); giving up on .LD.gz output\n";
        failed_ = true;
        return;
    }
    char name[64];
    std::snprintf(name, sizeof(name), "part%06d", part_seq_++);
    const std::string path = tmpdir_ + "/" + name;
    cur_ = gzopen(path.c_str(), "wb");
    if (!cur_) {
        std::cerr << "error: cannot open LD temp file " << path << std::endl;
        failed_ = true;
        return;
    }
    parts_.push_back(Part{label, path});
}

void LdWriter::close_part() {
    if (cur_) {
        if (gzclose(cur_) != Z_OK) failed_ = true;
        cur_ = nullptr;
    }
}

void LdWriter::emit_rr(const std::string& chr, int64_t p1, int64_t p2,
                       double rr, int64_t dis) {
    if (!cur_) return;
    char buf[128];
    if (!write_row(cur_, buf, sizeof(buf), "%s\t%lld\t%lld\t%.4f\t%lld\n",
                   chr.c_str(), static_cast<long long>(p1),
                   static_cast<long long>(p2), rr,
                   static_cast<long long>(dis)))
        failed_ = true;
}

void LdWriter::emit_d_rr(const std::string& chr, int64_t p1, int64_t p2,
                         double d, double rr, int64_t dis) {
    if (!cur_) return;
    char buf[160];
    if (!write_row(cur_, buf, sizeof(buf),
                   "%s\t%lld\t%lld\t%.4f\t%.4f\t%lld\n", chr.c_str(),
                   static_cast<long long>(p1), static_cast<long long>(p2), d, rr,
                   static_cast<long long>(dis)))
        failed_ = true;
}

void LdWriter::emit_d_lod_rr(const std::string& chr, int64_t p1, int64_t p2,
                             double d, double lod, double rr, int64_t dis) {
    if (!cur_) return;
    char buf[192];
    if (!write_row(cur_, buf, sizeof(buf),
                   "%s\t%lld\t%lld\t%.4f\t%.4f\t%.4f\t%lld\n", chr.c_str(),
                   static_cast<long long>(p1), static_cast<long long>(p2), d, lod,
                   rr, static_cast<long long>(dis)))
        failed_ = true;
}

void LdWriter::emit_d_lod_rr_ci(const std::string& chr, int64_t p1, int64_t p2,
                                double d, double lod, double rr, int low_i,
                                int high_i, int64_t dis) {
    if (!cur_) return;
    char buf[224];
    if (!write_row(cur_, buf, sizeof(buf),
                   "%s\t%lld\t%lld\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%lld\n",
                   chr.c_str(), static_cast<long long>(p1),
                   static_cast<long long>(p2), d, lod, rr, low_i / 100.0,
                   high_i / 100.0, static_cast<long long>(dis)))
        failed_ = true;
}

void LdWriter::absorb(LdWriter& other) {
    close_part();
    other.close_part();
    // A part that could not be written or closed in the worker is missing
    // rows; carrying that over keeps the failure from turning into a
    // silently short .LD.gz.
    if (other.failed_) failed_ = true;
    for (Part& p : other.parts_) {
        // Move the part file into our own temp dir so it outlives `other`
        // (which is destroyed when its worker goes out of scope).
        char dst_name[64];
        std::snprintf(dst_name, sizeof(dst_name), "part%06d", part_seq_++);
        const std::string dst = tmpdir_ + "/" + dst_name;
        bool ok = std::rename(p.path.c_str(), dst.c_str()) == 0;
        if (!ok) {
            // Cross-device: fall back to a byte copy.
            FILE* src = std::fopen(p.path.c_str(), "r");
            FILE* f = std::fopen(dst.c_str(), "w");
            ok = src && f;
            if (ok) {
                char buf[16384];
                size_t n;
                while ((n = std::fread(buf, 1, sizeof(buf), src)) > 0) {
                    if (std::fwrite(buf, 1, n, f) != n) {
                        ok = false;
                        break;
                    }
                }
            }
            if (src) std::fclose(src);
            if (f) std::fclose(f);
            if (!ok) {
                std::cerr << "warning: could not merge LD temp part (label "
                          << p.label << "); its rows are lost\n";
                failed_ = true;
                continue;
            }
        }
        parts_.push_back(Part{p.label, dst});
    }
    other.parts_.clear();
}

bool LdWriter::finish(const std::string& prefix) {
    close_part();
    const std::string out_path = prefix + ".LD.gz";
    if (failed_) {
        // A part was lost or a write failed: the table would be short, and a
        // short table that looks complete is worse than no table at all.
        std::cerr << "error: LD output incomplete (a part could not be written); "
                     "not writing " << out_path << std::endl;
        std::remove(out_path.c_str());
        return false;
    }
    FILE* out = std::fopen(out_path.c_str(), "wb");
    if (!out) {
        std::cerr << "open LD OUT File error: " << out_path << std::endl;
        return false;
    }

    std::sort(parts_.begin(), parts_.end(),
              [](const Part& x, const Part& y) { return x.label < y.label; });

    // Header as its own gzip member, then one member per part, in label order.
    // (Concatenated gzip members are a valid gzip stream whose decompressed
    // content is the concatenation, so the table content is the same as the
    // original's single-stream ogzstream output.)
    std::string hdr;
    if (!gzip_member(header_for(out_type_), hdr)) {
        std::cerr << "error: header gzip compression failed" << std::endl;
        std::fclose(out);
        std::remove(out_path.c_str());  // don't leave a truncated file behind
        return false;
    }
    if (std::fwrite(hdr.data(), 1, hdr.size(), out) != hdr.size()) {
        std::cerr << "write LD OUT File error: " << out_path << std::endl;
        std::fclose(out);
        std::remove(out_path.c_str());  // don't leave a truncated file behind
        return false;
    }

    char buf[16384];
    for (const Part& p : parts_) {
        FILE* f = std::fopen(p.path.c_str(), "rb");
        if (!f) {
            std::cerr << "error: cannot open LD temp part " << p.path << std::endl;
            std::fclose(out);
            std::remove(out_path.c_str());  // don't leave a truncated file behind
            return false;
        }
        size_t n;
        while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) {
            if (std::fwrite(buf, 1, n, out) != n) {
                std::cerr << "write LD OUT File error: " << out_path << std::endl;
                std::fclose(f);
                std::fclose(out);
                std::remove(out_path.c_str());  // don't leave a truncated file behind
                return false;
            }
        }
        std::fclose(f);
    }
    // fclose() flushes the last block; a failure here means the file on disk
    // is short, so it must not be published as a complete table.
    if (std::fclose(out) != 0) {
        std::cerr << "write LD OUT File error: " << out_path << std::endl;
        std::remove(out_path.c_str());  // don't leave a truncated file behind
        return false;
    }
    return true;
}

}  // namespace pld2
