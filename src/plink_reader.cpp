#include "plink_reader.h"
#include "textutil.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <thread>
#include <vector>

#include "region.h"

namespace pld2 {

namespace {

// PLINK .bed 2-bit genotype codes (variant-major .bed, magic 6c 1b 01):
// 00 homozygous A1, 01 missing, 10 heterozygous, 11 homozygous A2.
constexpr unsigned kHom1 = 0;
constexpr unsigned kMissing = 1;
constexpr unsigned kHet = 2;
constexpr unsigned kHom2 = 3;

// One .bim line. The allele strings are only needed for validation and for the
// -OutFilterSNP listing: the encoding is relative to the most frequent .bed
// code, exactly like the VCF path encodes relative to the most frequent GT
// character ('0' for a biallelic GT field).
struct BimRow {
    std::string chrom;
    std::string vid;
    std::string a1;
    std::string a2;
    int64_t pos = 0;
    bool indel = false;  // allele is not a single A/C/G/T base
    bool mono = false;   // a1 == a2
};

const std::vector<int>& empty_cols() {
    static const std::vector<int> kEmpty;
    return kEmpty;
}


bool is_acgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

// One uppercase single-base allele, or "" when the .bim allele is not a SNP
// base (multi-character indel, or PLINK's D/I/0 codes).
std::string normalized_allele(const std::string& a) {
    if (a.size() != 1) return std::string();
    const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(a[0])));
    return is_acgt(c) ? std::string(1, c) : std::string();
}

bool read_fam(const std::string& path, std::vector<std::string>& ids) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "open PLINK .fam IN File error: " << path << std::endl;
        return false;
    }
    std::set<std::string> seen;
    std::string line;
    std::vector<std::string> toks;  // reused across lines
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        split_ws_into(line, toks, /*cr=*/true);
        if (toks.size() < 2) continue;  // guard against malformed rows
        // Column 2 is the IID: plink --recode vcf names VCF samples after it.
        if (!seen.insert(toks[1]).second) {
            std::cerr << "Sample [ " << toks[1]
                      << " ] can be found [Repeat] in the PLINK .fam" << std::endl;
            return false;
        }
        ids.push_back(toks[1]);
    }
    if (ids.empty()) {
        std::cerr << "the PLINK .fam has no sample: " << path << std::endl;
        return false;
    }
    return true;
}

bool read_bim(const std::string& path, std::vector<BimRow>& rows) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "open PLINK .bim IN File error: " << path << std::endl;
        return false;
    }
    std::string line;
    std::vector<std::string> toks;  // reused across lines
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        split_ws_into(line, toks, /*cr=*/true);
        if (toks.size() < 6) continue;  // guard against malformed rows
        BimRow row;
        row.chrom = toks[0];
        row.vid = toks[1];
        row.pos = strtoll(toks[3].c_str(), nullptr, 10);
        row.a1 = normalized_allele(toks[4]);
        row.a2 = normalized_allele(toks[5]);
        // Not single-base alleles: same fate as a multi-character VCF REF/ALT
        // (Indel filter), and no .bed decoding needed.
        row.indel = row.a1.empty() || row.a2.empty();
        row.mono = !row.indel && row.a1 == row.a2;
        rows.push_back(std::move(row));
    }
    return true;
}

// Filters and encodes one variant, appending the packed site to every group's
// per-contig buffer. The group buffers are indexed by group, and (with -L) by
// the owning region label instead of the contig name. Returns true when the
// site was kept for the all-samples buffer (used by -OutFilterSNP).
bool process_variant(const BimRow& row, const unsigned char* bed,
                     const std::vector<SubPop>& subs, int nsamples,
                     const Options& opt, const RegionSet* rset,
                     std::vector<std::map<std::string, std::vector<GenoSite>>>& raw,
                     ReadStats& stats) {
    if (row.indel) {
        // Site-level, like the VCF Indel filter: counted once per variant.
        stats.skip_indel++;
        return false;
    }

    const std::string* key = &row.chrom;
    if (rset) {
        const std::string& label = rset->find(row.chrom, row.pos);
        if (label.empty()) return false;  // outside every -L region
        key = &label;
    }

    bool kept0 = false;
    for (size_t g = 0; g < raw.size(); ++g) {
        const std::vector<int>& cols = subs.empty() ? empty_cols() : subs[g].sample_cols;
        const int n_eff = subs.empty() ? nsamples : subs[g].n_samples;
        // With -SubPop the resolved columns index the synthetic
        // "#CHROM POS <IID>..." header built from the .fam, hence the -2.
        auto sample_index = [&](int k) -> size_t {
            return cols.empty() ? static_cast<size_t>(k)
                                : static_cast<size_t>(cols[k]) - 2;
        };

        // Allele counts over the pseudo GT characters: .bed code kHom1 is '0',
        // kHom2 is '1', kHet contributes one of each. This is the PLINK
        // counterpart of counting the raw GT characters in the VCF path, and it
        // makes the best/sed selection ('0' wins ties) come out identical.
        int64_t n0 = 0;
        int64_t n1 = 0;
        int het_count = 0;
        int miss_count = 0;
        for (int k = 0; k < n_eff; ++k) {
            const size_t idx = sample_index(k);
            const unsigned code = (bed[idx >> 2] >> (2 * (idx & 3))) & 3u;
            if (code == kMissing) {
                ++miss_count;
            } else if (code == kHet) {
                ++het_count;
                ++n0;
                ++n1;
            } else if (code == kHom1) {
                n0 += 2;
            } else {
                n1 += 2;
            }
        }

        // Filters, in the VCF reader's order and with its formulas.
        if ((miss_count * 1.0 / n_eff) > opt.miss) continue;
        if ((het_count * 1.0 / n_eff) > opt.het) continue;
        const int base_count = (n0 > 0 ? 1 : 0) + (n1 > 0 ? 1 : 0);
        if (row.mono || base_count == 1 || base_count > 2) {
            stats.skip_non_biallelic++;
            continue;
        }
        const int64_t max_c = n1 > n0 ? n1 : n0;
        const int64_t sed_c = n1 > n0 ? n0 : n1;
        if ((sed_c * 1.0 / (sed_c + max_c)) < opt.maf) continue;

        // Encode 0 = best, 1 = second, 2 = missing. A heterozygote is the
        // (best, sed) pseudo-haplotype pair, i.e. always (0, 1).
        const uint8_t hom1_code = (n1 > n0) ? 1 : 0;
        const uint8_t hom2_code = static_cast<uint8_t>(1 - hom1_code);
        std::vector<uint8_t> vals;
        vals.reserve(2 * static_cast<size_t>(n_eff));
        for (int k = 0; k < n_eff; ++k) {
            const size_t idx = sample_index(k);
            const unsigned code = (bed[idx >> 2] >> (2 * (idx & 3))) & 3u;
            uint8_t v0, v1;
            if (code == kHet) {
                v0 = 0;
                v1 = 1;
            } else if (code == kHom1) {
                v0 = v1 = hom1_code;
            } else if (code == kHom2) {
                v0 = v1 = hom2_code;
            } else {
                v0 = v1 = 2;
            }
            vals.push_back(v0);
            vals.push_back(v1);
        }

        std::vector<GenoSite>& sites = raw[g][*key];
        sites.push_back(GenoSite{row.pos, std::move(vals)});
        stats.sites_kept++;
        if (g == 0) kept0 = true;
    }
    return kept0;
}

}  // namespace

int read_plink(const std::string& path, const Options& opt,
               const std::string& prefix,
               const std::vector<SubGroupSpec>& specs, ReadStats& stats,
               std::vector<std::vector<GenoBuf>>& out) {
    std::string base = path;
    if (base.size() > 4 && base.compare(base.size() - 4, 4, ".bed") == 0) {
        base.resize(base.size() - 4);
    }
    const std::string bed_path = base + ".bed";
    const std::string bim_path = base + ".bim";
    const std::string fam_path = base + ".fam";

    std::vector<std::string> ids;
    if (!read_fam(fam_path, ids)) return 0;
    const int nsamples = static_cast<int>(ids.size());

    std::vector<BimRow> rows;
    if (!read_bim(bim_path, rows)) return 0;
    const size_t nbytes = (static_cast<size_t>(nsamples) + 3) / 4;

    FILE* fp = fopen(bed_path.c_str(), "rb");
    if (!fp) {
        std::cerr << "open PLINK .bed IN File error: " << bed_path << std::endl;
        return 0;
    }
    // The .bed magic is THREE bytes, not two: 6c 1b/1a selects the layout and
    // the third byte selects the mode. 6c 1b 01 is the only variant-major,
    // SNP-only, 2-bit layout this reader supports; 6c 1b 00 stores SNPs as 2-bit
    // and indels as two 2-bit hardcalls and must not be misread as pure 2-bit.
    unsigned char magic[3] = {0, 0, 0};
    if (fread(magic, 1, 3, fp) != 3) {
        std::cerr << "not a PLINK .bed file: " << bed_path << " (too short)" << std::endl;
        fclose(fp);
        return 0;
    }
    if (magic[0] == 0x6c && magic[1] == 0x1a) {
        std::cerr << "the PLINK .bed is individual-major (6c 1a); regenerate the "
                     "dataset with plink --make-bed"
                  << std::endl;
        fclose(fp);
        return 0;
    }
    if (magic[0] != 0x6c || magic[1] != 0x1b) {
        std::cerr << "not a PLINK .bed file: " << bed_path << " (bad magic bytes)"
                  << std::endl;
        fclose(fp);
        return 0;
    }
    if (magic[2] != 0x01) {
        std::cerr << "the PLINK .bed mixes SNPs and indels (magic 6c 1b "
                  << static_cast<int>(magic[2])
                  << "); this reader supports the SNP-only layout (6c 1b 01) "
                     "only, so regenerate the dataset without indels"
                  << std::endl;
        fclose(fp);
        return 0;
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        std::cerr << "the PLINK .bed is not seekable: " << bed_path << std::endl;
        fclose(fp);
        return 0;
    }
    const long file_size = ftell(fp);
    // The 3-byte magic sits in front of the genotype plane, so the expected
    // length is 3 + nvariants * ceil(nsamples/4).
    const long long expect_size =
        3 + static_cast<long long>(rows.size()) * static_cast<long long>(nbytes);
    if (file_size != expect_size) {
        // The usual cause is PLINK's >524288-variants-per-chromosome layout,
        // which stores one full plane per variant instead of 2-bit codes.
        std::cerr << "the PLINK .bed size does not match the .bim: " << file_size
                  << " bytes, expected " << expect_size << " for " << rows.size()
                  << " variants x " << nsamples << " samples; regenerate the "
                     "dataset with plink --make-bed"
                  << std::endl;
        fclose(fp);
        return 0;
    }
    fclose(fp);

    // -SubPop: the synthetic "#CHROM POS <IID>..." header makes the shared
    // resolver usable; sample columns then index the .fam at col - 2.
    const size_t nbufs = specs.empty() ? 1 : specs.size();
    std::vector<SubPop> subs;
    if (!specs.empty()) {
        std::vector<std::string> header;
        header.push_back("#CHROM");
        header.push_back("POS");
        header.insert(header.end(), ids.begin(), ids.end());
        if (!resolve_subpop_header(header, specs, "[found in PLINK .fam]",
                                   "in the PLINK .fam", 2, subs)) {
            return 0;
        }
    }

    // -L: keep the sites of each region in its own buffer, so pairs across
    // regions never meet (identical to the read_vcf fallback).
    RegionSet rset(opt.regions);
    const RegionSet* rsetp = opt.regions.empty() ? nullptr : &rset;

    // -OutFilterSNP lists the passing variants in .bim order, so it keeps the
    // single-threaded pass.
    FILE* filt = nullptr;
    if (opt.out_filter_snp && specs.empty()) {
        const std::string filt_path = prefix + ".plink.filter.bim";
        filt = fopen(filt_path.c_str(), "w");
        if (!filt) {
            std::cerr << "open OUT File error: " << filt_path << std::endl;
            return 0;
        }
    }

    // Decode in -T chunks over the variant range: the .bed is random access, so
    // the only shared state is the read-only .bim.
    const size_t nrows = rows.size();
    size_t nworkers = 1;
    if (!filt && opt.threads > 1 && nrows >= 4096) {
        const size_t kRowsPerWorker = 4096;
        nworkers = std::min<size_t>(static_cast<size_t>(opt.threads),
                                    (nrows + kRowsPerWorker - 1) / kRowsPerWorker);
        nworkers = std::min<size_t>(nworkers, 64);
    }

    struct Chunk {
        std::vector<std::map<std::string, std::vector<GenoSite>>> raw;
        ReadStats stats;
        bool failed = false;
    };
    std::vector<Chunk> chunks(nworkers);
    for (Chunk& c : chunks) c.raw.resize(nbufs);

    auto run_chunk = [&](size_t w, size_t begin, size_t end) {
        FILE* f = fopen(bed_path.c_str(), "rb");
        if (!f) {
            std::cerr << "open PLINK .bed IN File error: " << bed_path << std::endl;
            chunks[w].failed = true;
            return;
        }
        setvbuf(f, nullptr, _IOFBF, 1 << 20);
        // Skip the 3-byte magic; genotype data for variant `begin` starts at
        // 3 + begin * nbytes. Using 2 here shifted every chunk one byte early
        // and silently decoded the wrong genotypes.
        if (fseek(f, static_cast<long>(3 + begin * nbytes), SEEK_SET) != 0) {
            fclose(f);
            chunks[w].failed = true;
            return;
        }
        std::vector<unsigned char> buf(nbytes);
        for (size_t i = begin; i < end; ++i) {
            if (fread(buf.data(), 1, nbytes, f) != nbytes) {
                std::cerr << "read PLINK .bed error: " << bed_path << std::endl;
                fclose(f);
                chunks[w].failed = true;
                return;
            }
            chunks[w].stats.sites_in++;
            if (process_variant(rows[i], buf.data(), subs, nsamples, opt, rsetp,
                                chunks[w].raw, chunks[w].stats) &&
                filt) {
                fprintf(filt, "%s\t%s\t0\t%lld\t%s\t%s\n", rows[i].chrom.c_str(),
                        rows[i].vid.c_str(), static_cast<long long>(rows[i].pos),
                        rows[i].a1.c_str(), rows[i].a2.c_str());
            }
        }
        fclose(f);
    };

    std::vector<std::thread> pool;
    pool.reserve(nworkers - 1);
    for (size_t w = 1; w < nworkers; ++w) {
        const size_t begin = nrows * w / nworkers;
        const size_t end = nrows * (w + 1) / nworkers;
        pool.emplace_back(run_chunk, w, begin, end);
    }
    run_chunk(0, 0, nrows / nworkers);
    for (std::thread& t : pool) t.join();

    if (filt) fclose(filt);
    for (const Chunk& c : chunks) {
        if (c.failed) return 0;
    }

    for (const Chunk& c : chunks) {
        stats.sites_in += c.stats.sites_in;
        stats.sites_kept += c.stats.sites_kept;
        stats.skip_indel += c.stats.skip_indel;
        stats.skip_non_biallelic += c.stats.skip_non_biallelic;
    }

    // Concatenate the chunks in worker order, so each contig's site list is in
    // .bim order and merge_geno_buffers' stable sort/dedupe keeps the same
    // first occurrence as the single-threaded pass.
    out.assign(nbufs, std::vector<GenoBuf>());
    for (size_t g = 0; g < nbufs; ++g) {
        std::map<std::string, std::vector<GenoSite>> merged;
        for (size_t w = 0; w < nworkers; ++w) {
            for (auto& kv : chunks[w].raw[g]) {
                std::vector<GenoSite>& dst = merged[kv.first];
                dst.insert(dst.end(), std::make_move_iterator(kv.second.begin()),
                           std::make_move_iterator(kv.second.end()));
            }
        }
        std::vector<GenoBuf> gbs;
        gbs.reserve(merged.size());
        for (auto& kv : merged) gbs.push_back(GenoBuf{kv.first, std::move(kv.second)});
        out[g] = std::move(gbs);
    }
    return 1;
}

}  // namespace pld2