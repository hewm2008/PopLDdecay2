#include "genotype_reader.h"

#include <zlib.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>  // hts_get_bgzfp: the threaded decompression gate

#include "subpop.h"
#include "textutil.h"

namespace pld2 {

namespace {

// Tokenization on space/tab (replicates split(line,inf," \t")).

// IUPAC code -> allele pair, as in FilterGenotype.h / FileDeal.h (unknown
// tokens map to "").
// The IUPAC code table: the original looked each token up in a
// map<string,string>. Every key there is a single character, so a 256-entry
// table indexed by the first character answers the same question -- and a
// token of any other length is in no entry either, so it maps to "" (and
// counts as a missing call).
struct IupacTable {
    const char* code[256];
    IupacTable() {
        for (int i = 0; i < 256; ++i) code[i] = "";
        code[(unsigned char)'M'] = "AC";
        code[(unsigned char)'K'] = "GT";
        code[(unsigned char)'Y'] = "CT";
        code[(unsigned char)'R'] = "AG";
        code[(unsigned char)'W'] = "AT";
        code[(unsigned char)'S'] = "CG";
        code[(unsigned char)'C'] = "CC";
        code[(unsigned char)'G'] = "GG";
        code[(unsigned char)'T'] = "TT";
        code[(unsigned char)'A'] = "AA";
        code[(unsigned char)'-'] = "NN";
        code[(unsigned char)'N'] = "NN";
    }
};
const IupacTable kIupac;

// The allele pair a token stands for ("" when the token is unknown).
inline const char* iupac_code(const std::string& tok) {
    return tok.size() == 1 ? kIupac.code[(unsigned char)tok[0]] : "";
}

// Allele pair -> IUPAC code, as in FileDeal.h Read_Genotype_IN.
const std::map<std::string, char>& fwd_map() {
    static const std::map<std::string, char> m = {
        {"AC", 'M'}, {"CA", 'M'}, {"GT", 'K'}, {"TG", 'K'}, {"CT", 'Y'},
        {"TC", 'Y'}, {"AG", 'R'}, {"GA", 'R'}, {"AT", 'W'}, {"TA", 'W'},
        {"CG", 'S'}, {"GC", 'S'}, {"AA", 'A'}, {"TT", 'T'}, {"CC", 'C'},
        {"GG", 'G'}};
    return m;
}

// Allele counts for one row: a 256-entry array plus the list of characters
// that actually occurred. This replaces the original's map<char,int> -- the
// counts are the same, and the characters are walked in ascending order, so
// the best/second-best selection is identical -- but it costs no allocation
// per row.
struct AlleleCounts {
    int count[256] = {0};
    unsigned char touched[256] = {0};
    int ntouched = 0;

    void add(char c) {
        const unsigned char u = static_cast<unsigned char>(c);
        if (count[u]++ == 0) touched[ntouched++] = u;
    }
    void reset() {
        for (int i = 0; i < ntouched; ++i) count[touched[i]] = 0;
        ntouched = 0;
    }
    // Ascending character order, like map<char,int>'s iteration.
    void sort_touched() {
        for (int i = 1; i < ntouched; ++i) {
            const unsigned char v = touched[i];
            int j = i - 1;
            while (j >= 0 && touched[j] > v) {
                touched[j + 1] = touched[j];
                --j;
            }
            touched[j + 1] = v;
        }
    }
};

// best / second-best allele selection over a char-count map, replicating the
// original `map<char,int>::begin()` loops ('N' excluded; later chars win ties
// for the second-best). `miss` receives the 'N' count.
struct Best {
    char best = 'N';
    char sed = 'N';
    int best_c = 0;
    int sed_c = 0;
    int base_count = 0;
    int miss = 0;
};

Best pick_best(const AlleleCounts& Count) {
    Best r;
    for (int i = 0; i < Count.ntouched; ++i) {
        const unsigned char c = Count.touched[i];
        const int n = Count.count[c];
        if (c == 'N') {
            r.miss = n;
            continue;
        }
        if (n > r.best_c) {
            r.sed_c = r.best_c;
            r.sed = r.best;
            r.best_c = n;
            r.best = static_cast<char>(c);
        } else if (n >= r.sed_c) {
            r.sed_c = n;
            r.sed = static_cast<char>(c);
        }
        r.base_count++;
    }
    return r;
}

// Per-token counting replicating FilterGeno: unknown tokens look up "" and
// count as '\0' (the libstdc++ SSO read of the null terminator).
void count_tokens(const std::vector<std::string>& toks, int* het,
                  AlleleCounts& Count) {
    *het = 0;
    for (size_t ii = 2; ii < toks.size(); ++ii) {
        const char* A_tmp = iupac_code(toks[ii]);
        char c0 = A_tmp[0] ? A_tmp[0] : '\0';
        char c1 = A_tmp[0] && A_tmp[1] ? A_tmp[1] : '\0';
        Count.add(c0);
        Count.add(c1);
        if (c0 != c1) (*het)++;
    }
    Count.sort_touched();
}

// FilterGeno filter for one row (Cut3base is always on in the main tool).
// Fills `b` with the allele counts, so the caller does not count the row
// a second time for the encoding.
bool filter_geno(const std::vector<std::string>& toks, const Options& opt,
                 Best& b) {
    int het = 0;
    AlleleCounts Count;
    count_tokens(toks, &het, Count);
    b = pick_best(Count);
    const int sample = static_cast<int>(toks.size()) - 2;
    if ((b.miss * 1.0) / (sample * 2.0) > opt.miss) return false;
    if ((het * 1.0) / sample > opt.het) return false;
    if ((b.sed_c * 1.0) / (sample * 2.0) < opt.maf) return false;
    if (b.base_count < 2) return false;
    if (b.base_count > 2) return false;
    return true;
}

// GetBestBase encoding (replicates FileDeal.h GetBestBase: allele pair -> IUPAC
// lookup for the het code, per-token Allele2double resolution, unknown tokens
// coerced to Het_base -> (0,1)).
std::vector<uint8_t> encode_geno(const std::vector<std::string>& toks,
                                 const Best& b) {
    const std::string A_base(1, b.best);
    const std::string B_base(1, b.sed);
    const std::string Het_base = A_base + B_base;
    auto it_sa = fwd_map().find(Het_base);
    const std::string C_base =
        (it_sa != fwd_map().end()) ? std::string(1, it_sa->second) : "N";

    // Allele2double: the original built a map per row with the entries
    // "-"->"NN", "N"->"NN", "n"->"NN", A_base, B_base and C_base, and a
    // token that is in no entry fell back to Het_base. A later entry wins if
    // two keys are the same, hence the order below.
    const std::string AA = A_base + A_base;
    const std::string BB = B_base + B_base;
    std::vector<uint8_t> vals;
    vals.reserve(2 * (toks.size() - 2));
    for (size_t ii = 2; ii < toks.size(); ++ii) {
        const std::string& tok = toks[ii];
        std::string A_tmp;
        if (tok == C_base) {
            A_tmp = Het_base;
        } else if (tok == B_base) {
            A_tmp = BB;
        } else if (tok == A_base) {
            A_tmp = AA;
        } else if (tok == "n" || tok == "N" || tok == "-") {
            A_tmp = "NN";
        } else {
            A_tmp = Het_base;  // unknown token
        }
        const char c0 = A_tmp.empty() ? '\0' : A_tmp[0];
        const char c1 = A_tmp.size() > 1 ? A_tmp[1] : '\0';
        vals.push_back((uint8_t)((c0 == b.best) ? 0 : (c0 == b.sed) ? 1 : 2));
        vals.push_back((uint8_t)((c1 == b.best) ? 0 : (c1 == b.sed) ? 1 : 2));
    }
    return vals;
}

// Decodes one packed site back to per-haplotype values (for cross-input merges).
std::vector<uint8_t> decode_site(const SnpData& sd, size_t n_haps) {
    std::vector<uint8_t> vals(n_haps, 0);
    for (size_t h = 0; h < n_haps; ++h) {
        const size_t w = h / 64;
        const size_t bit = h % 64;
        if (w < sd.miss.size() && (sd.miss[w] & (uint64_t(1) << bit)))
            vals[h] = 2;
        else if (w < sd.minor.size() && (sd.minor[w] & (uint64_t(1) << bit)))
            vals[h] = 1;
    }
    return vals;
}

// Packs a site's values into bit-packed words at the contig's n_haps/nwords
// (shorter sites padded with best 0, longer truncated). Words are slices of
// `arena` (owned by the ContigData being built).
SnpData pack_site(const GenoSite& gs, size_t n_haps, size_t nwords,
                  SnpArena& arena) {
    SnpData sd;
    sd.pos = gs.pos;
    sd.minor.fill(arena, nwords, 0);
    sd.miss.fill(arena, nwords, 0);
    const size_t m = std::min(n_haps, gs.vals.size());
    for (size_t h = 0; h < m; ++h) {
        const uint8_t v = gs.vals[h];
        const size_t w = h / 64;
        const size_t bit = h % 64;
        if (v == 1)
            sd.minor[w] |= (uint64_t(1) << bit);
        else if (v == 2)
            sd.miss[w] |= (uint64_t(1) << bit);
    }
    site_counts(sd.minor.data(), sd.miss.data(), nwords, sd.minor_count,
                sd.has_missing);
    return sd;
}

// ---- non-SubPop path: FilterGeno + GetBestBase + optional filter.gz ----
int read_genotype_full(const std::string& path, const Options& opt,
                       const std::string& prefix, std::vector<GenoBuf>& out,
                       ReadStats& stats) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "open Genotype File IN File error: " << path << std::endl;
        return 0;
    }
    if (opt.threads > 1 && hts_get_bgzfp(fp)) hts_set_threads(fp, opt.threads);

    // OutFilterSNP output: the passing data lines verbatim (plain gzip,
    // byte-compatible with the original ogzstream FilterGeno output).
    gzFile filt = nullptr;
    std::string filt_path;
    if (opt.out_filter_snp) {
        filt_path = prefix + ".genotype.filter.gz";
        filt = gzopen(filt_path.c_str(), "wb");
        if (!filt) {
            std::cerr << "open OUT File error: " << filt_path << std::endl;
            hts_close(fp);
            return 0;
        }
    }

    std::map<std::string, std::vector<GenoSite>> raw;
    kstring_t line = {0, 0, 0};
    std::vector<std::string> toks;  // reused across lines
    int ret;
    while ((ret = hts_getline(fp, '\n', &line)) >= 0) {
        if (line.l == 0) continue;
        if (line.s[0] == '#') continue;
        split_ws_into(line.s, line.l, toks);
        if (toks.size() < 3) continue;  // guard against malformed rows
        stats.sites_in++;
        Best b;
        if (!filter_geno(toks, opt, b)) continue;
        stats.sites_kept++;
        if (filt) {
            gzwrite(filt, line.s, static_cast<unsigned>(line.l));
            gzputc(filt, '\n');
        }
        std::vector<uint8_t> vals = encode_geno(toks, b);
        const int64_t pos = strtoll(toks[1].c_str(), nullptr, 10);
        raw[toks[0]].push_back(GenoSite{pos, std::move(vals)});
    }
    if (filt) gzclose(filt);
    hts_close(fp);

    for (auto& [name, sites] : raw) {
        out.push_back(GenoBuf{name, std::move(sites)});
    }
    return 1;
}

// ---- SubPop path: Read_SubPopGenotype_IN inline filter + encode ----
int read_genotype_subpop(const std::string& path, const Options& opt,
                         const std::vector<SubGroupSpec>& specs,
                         std::vector<std::vector<GenoBuf>>& out,
                         ReadStats& stats) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "open Genotype File IN File error: " << path << std::endl;
        return 0;
    }

    // header -> matched sample columns (replicates the #CHROM prelude,
    // including the "wrong Line" abort when data precedes the header)
    std::vector<std::string> header;
    kstring_t line = {0, 0, 0};
    int ret;
    while ((ret = hts_getline(fp, '\n', &line)) >= 0) {
        if (line.l == 0) continue;
        std::string s(line.s, line.l);
        if (s[0] == '#' && s.size() > 1 && s[1] == '#') continue;
        if (s[0] == '#') {
            if (s.size() > 1 && s[1] != '#') {
                std::vector<std::string> inf = split_ws(s);
                if (inf.empty() || inf[0] != "#CHROM") continue;
                header = inf;
                break;
            }
            continue;  // '#' alone: old UB path, treat as non-header comment
        }
        if (s.size() > 1 && s[1] == '#') continue;  // old: falls through all branches
        std::cerr << "wrong Line : " << s << std::endl;
        std::cerr << "Genotype Header something wrong, can't find sample info before site info"
                  << std::endl;
        std::cerr << "iTools   Formtools  VCF2Genotype  -InPut  in.vcf  -OutPut  out.genotype  -WithHeader   -NoRef"
                  << std::endl;
        std::cerr << "Genotype Header sample info Flag : [  #CHROM  ] " << std::endl;
        hts_close(fp);
        return 0;
    }

    std::vector<SubPop> subs;
    if (!resolve_subpop_header(header, specs, "[found in Genotype]",
                               "in the Genotype Header", 2, subs)) {
        hts_close(fp);
        return 0;
    }
    const size_t nbufs = subs.size();

    std::vector<std::map<std::string, std::vector<GenoSite>>> raws(nbufs);
    int bad_site = 0;
    std::vector<std::string> inf;  // reused across lines
    while ((ret = hts_getline(fp, '\n', &line)) >= 0) {
        if (line.l == 0) continue;
        // Only empty lines are skipped here, like the original's data loop
        // (a '#' line after the header is not special there).
        split_ws_into(line.s, line.l, inf);
        if (inf.size() < 2) continue;  // guard against malformed rows
        stats.sites_in++;
        bool kept_by_any = false;

        for (size_t g = 0; g < nbufs; ++g) {
            const int number_sub = subs[g].n_samples;
            const std::vector<int>& sample_site = subs[g].sample_cols;

            AlleleCounts Count;
            int het_c = 0;
            int miss_c = 0;
            for (int kk = 0; kk < number_sub; ++kk) {
                const size_t idx = static_cast<size_t>(sample_site[kk]);
                const char* G = iupac_code(idx < inf.size() ? inf[idx] : "");
                const char ABase = G[0] ? G[0] : '\0';
                if (ABase == 'N') {
                    miss_c++;
                } else {
                    const char BBase = (G[0] && G[1]) ? G[1] : '\0';
                    if (ABase != BBase) het_c++;
                    Count.add(ABase);
                    Count.add(BBase);
                }
            }
            Count.sort_touched();

            if ((miss_c * 1.0 / number_sub) > opt.miss) continue;
            if ((het_c * 1.0 / number_sub) > opt.het) continue;

            Best b = pick_best(Count);
            if (b.base_count == 1 || b.base_count > 2) {
                bad_site++;
                continue;
            }
            if ((b.sed_c * 1.0 / (b.sed_c + b.best_c)) < opt.maf) continue;

            std::vector<uint8_t> vals;
            vals.reserve(2 * static_cast<size_t>(number_sub));
            for (int kk = 0; kk < number_sub; ++kk) {
                const size_t idx = static_cast<size_t>(sample_site[kk]);
                const char* G = iupac_code(idx < inf.size() ? inf[idx] : "");
                const char ABase = G[0] ? G[0] : '\0';
                char v0, v1;
                if (ABase == 'N') {
                    v0 = 'N';
                    v1 = 'N';
                } else {
                    const char BBase = (G[0] && G[1]) ? G[1] : '\0';
                    if (ABase != BBase) {
                        v0 = b.best;
                        v1 = b.sed;
                    } else {
                        v0 = ABase;
                        v1 = BBase;
                    }
                }
                vals.push_back((uint8_t)((v0 == b.best) ? 0 : (v0 == b.sed) ? 1 : 2));
                vals.push_back((uint8_t)((v1 == b.best) ? 0 : (v1 == b.sed) ? 1 : 2));
            }
            const int64_t pos = strtoll(inf[1].c_str(), nullptr, 10);
            raws[g][inf[0]].push_back(GenoSite{pos, std::move(vals)});
            kept_by_any = true;
        }
        if (kept_by_any) stats.sites_kept++;
    }
    hts_close(fp);

    if (bad_site != 0) {
        std::cout << "Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"
                  << bad_site << std::endl;
    }

    out.clear();
    out.reserve(nbufs);
    for (size_t g = 0; g < nbufs; ++g) {
        std::vector<GenoBuf> gb;
        for (auto& [name, sites] : raws[g]) {
            gb.push_back(GenoBuf{name, std::move(sites)});
        }
        out.push_back(std::move(gb));
    }
    return 1;
}

}  // namespace

int read_genotype(const std::string& path, const Options& opt,
                  const std::string& prefix,
                  const std::vector<SubGroupSpec>& specs,
                  std::vector<std::vector<GenoBuf>>& out, ReadStats& stats) {
    if (specs.empty()) {
        out.clear();
        out.resize(1);
        return read_genotype_full(path, opt, prefix, out[0], stats);
    }
    return read_genotype_subpop(path, opt, specs, out, stats);
}

void merge_geno_buffers(std::vector<GenoBuf>& gbufs, std::vector<ContigData>& out) {
    if (gbufs.empty()) return;

    // name -> combined site list (existing contigs first, then genotype sites)
    std::map<std::string, std::vector<GenoSite>> merged;
    for (const ContigData& cd : out) {
        for (const SnpData& sd : cd.snps) {
            merged[cd.name].push_back(GenoSite{sd.pos, decode_site(sd, cd.n_haplotypes)});
        }
    }
    for (const GenoBuf& gb : gbufs) {
        std::vector<GenoSite>& v = merged[gb.name];
        v.insert(v.end(), gb.sites.begin(), gb.sites.end());
    }

    out.clear();
    out.reserve(merged.size());
    for (auto& [name, sites] : merged) {
        std::stable_sort(sites.begin(), sites.end(),
                         [](const GenoSite& a, const GenoSite& b) {
                             return a.pos < b.pos;
                         });
        std::vector<GenoSite> uniq;
        uniq.reserve(sites.size());
        for (auto& gs : sites) {
            if (!uniq.empty() && uniq.back().pos == gs.pos) continue;  // keep first
            uniq.push_back(std::move(gs));
        }
        if (uniq.empty()) continue;
        const size_t n_haps = uniq.front().vals.size();
        const size_t nwords = (n_haps + 63) / 64;
        ContigData cd;
        cd.name = name;
        cd.nwords = nwords;
        cd.n_haplotypes = n_haps;
        cd.snps.reserve(uniq.size());
        for (auto& gs : uniq) cd.snps.push_back(pack_site(gs, n_haps, nwords, cd.arena));
        out.push_back(std::move(cd));
    }
}

}  // namespace pld2