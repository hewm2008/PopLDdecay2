#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "region.h"

namespace pld2 {

// Per-contig allocation arena for the bit-packed SNP words. Slices are handed
// out from chunked storage, so a returned pointer never moves (chunks are
// never reallocated/moved). Owned by the ContigData that built the sites; the
// owning ContigData is move-only, so slices stay valid across vector moves.
struct SnpArena {
    static constexpr size_t kChunkWords = 1 << 16;  // 512 KiB of uint64_t
    std::vector<std::unique_ptr<uint64_t[]>> chunks;
    size_t chunk_used = 0;  // words used in the current chunk

    uint64_t* alloc(size_t nwords) {
        if (nwords == 0) return nullptr;
        if (nwords > kChunkWords) {
            // Oversized slice (only reachable for >~2.1M samples): give it its
            // own dedicated chunk sized to fit, so it never overflows a normal
            // chunk. After it the next alloc() starts a fresh normal chunk.
            chunks.push_back(std::make_unique<uint64_t[]>(nwords));
            chunk_used = 0;
            return chunks.back().get();
        }
        if (chunks.empty() || chunk_used + nwords > kChunkWords) {
            chunks.push_back(std::make_unique<uint64_t[]>(kChunkWords));
            chunk_used = 0;
        }
        uint64_t* p = chunks.back().get() + chunk_used;
        chunk_used += nwords;
        return p;
    }
};

// Bit-packed words for one SNP (one bit per haplotype). The single-word case
// (n==1, <=64 haplotypes / <=32 samples) is stored inline in `one`; n>1 words
// are a slice of the owning ContigData's SnpArena (never heap-allocated per
// SNP). The struct is exactly 24 bytes (same as a std::vector) so SnpData
// keeps its previous size. `p` is a borrowed pointer -- the arena owns the
// storage.
struct SnpBits {
    size_t n = 0;
    uint64_t one = 0;
    uint64_t* p = nullptr;  // n==1 ? &one : arena slice

    size_t size() const { return n; }
    const uint64_t* data() const { return n == 1 ? &one : p; }
    uint64_t operator[](size_t i) const { return n == 1 ? one : p[i]; }
    uint64_t& operator[](size_t i) { return n == 1 ? one : p[i]; }

    // Copies `src` into an arena slice (or inline `one` for n==1).
    void assign(SnpArena& a, const uint64_t* src, size_t count) {
        n = count;
        if (count == 1) {
            one = src[0];
            p = nullptr;
        } else if (count > 1) {
            p = a.alloc(count);
            std::memcpy(p, src, count * sizeof(uint64_t));
        } else {
            p = nullptr;
        }
    }
    // Fills a fresh arena slice with `val`.
    void fill(SnpArena& a, size_t count, uint64_t val) {
        n = count;
        if (count == 1) {
            one = val;
            p = nullptr;
        } else if (count > 1) {
            p = a.alloc(count);
            std::fill(p, p + count, val);
        } else {
            p = nullptr;
        }
    }
};

// One SNP on a contig, bit-packed haplotypes.
// Per haplotype 2-bit semantics (replicates original Read_VCF_IN encoding):
//   0 = major (best), 1 = minor (second), 2 = missing
// minor[] bit set  => allele == 1
// miss[]  bit set  => allele == 2
struct SnpData {
    int64_t pos = 0;
    SnpBits minor;   // nwords words, one bit per haplotype
    SnpBits miss;    // nwords words, one bit per haplotype

    // SnpBits stores a borrowed arena pointer: a copy would alias the same
    // arena slices, so copies are rejected (moves only). SnpArena's
    // unique_ptr members already make ContigData move-only.
    SnpData() = default;
    SnpData(const SnpData&) = delete;
    SnpData& operator=(const SnpData&) = delete;
    SnpData(SnpData&&) = default;
    SnpData& operator=(SnpData&&) = default;
};

struct ContigData {
    std::string name;
    std::vector<SnpData> snps;     // sorted by pos (ascending), positions unique
    SnpArena arena;                // owns the bit-packed SNP words (see above)
    size_t nwords = 0;             // ceil(2*N / 64)
    size_t n_haplotypes = 0;       // 2*N
};

// Distance bins, indexed by physical distance in bp (1 .. max_dist_bp).
struct Bin {
    int64_t count = 0;
    double sumRR = 0.0;
    double sumD = 0.0;
};

struct Options {
    std::string in_vcf;
    std::string in_genotype;   // M6: -InGenotype input (native IUPAC format)
    std::string out_stat;
    std::vector<std::string> subpop_files; // M3c: -SubPop sample lists (repeatable)
    std::string ehh;           // M5: -EHH input
    std::string region;        // M3: raw -L argument
    std::vector<Region> regions;  // M3: parsed -L regions (sorted by label)
    int out_type = 1;          // OutType 0-8
    int max_dist_kb = 300;
    double maf = 0.005;
    double het = 0.88;
    double miss = 0.25;
    int method = 1;            // only Method 1 supported
    bool out_filter_snp = false; // M3: -OutFilterSNP
    int threads = 4;           // M2: parallel LD workers (>=1)
};

struct ReadStats {
    long long sites_in = 0;    // total records seen
    long long sites_kept = 0;  // records after filters
    long long skip_indel = 0;
    long long skip_non_biallelic = 0;
};

}  // namespace pld2