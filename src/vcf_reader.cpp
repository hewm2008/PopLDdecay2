#include "vcf_reader.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <zlib.h>

#include "cli.h"
#include "ldkernel.h"
#include "subpop.h"

namespace pld2 {

namespace {

// Header tokenization on space/tab (replicates split(line,inf," \t")).
std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> toks;
    std::string cur;
    for (char c : line) {
        if (c == ' ' || c == '\t') {
            if (!cur.empty()) {
                toks.push_back(cur);
                cur.clear();
            }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) toks.push_back(cur);
    return toks;
}

// Reads raw allele chars at index 0 and 2 of a token (replicates spline[jj][0]
// and spline[jj][2] reads, including the '0/.' / '.' counting quirks).
// Tab offsets are computed into a caller-owned scratch vector (reused across
// lines -> no per-line heap allocation). field() is only used where a real
// std::string is required (the chromosome map key); pos/ref/alt are read via
// data()/length() to avoid substring allocations.
class LineTokens {
  public:
    LineTokens(const char* s, size_t len, std::vector<size_t>& offs)
        : s_(s), len_(len), offs_(offs) {
        offs_.clear();
        size_t start = 0;
        for (size_t i = 0; i < len_; ++i) {
            if (s_[i] == '\t') {
                offs_.push_back(start);
                start = i + 1;
            }
        }
        offs_.push_back(start);
    }

    size_t off(int j) const {
        return (j < static_cast<int>(offs_.size())) ? offs_[j] : len_;
    }

    size_t end(int j) const {
        size_t o = off(j);
        if (j + 1 < static_cast<int>(offs_.size())) return offs_[j + 1] - 1;
        size_t e = len_;
        if (e > 0 && s_[e - 1] == '\n') e--;
        if (e > 0 && s_[e - 1] == '\r') e--;
        return (e > o) ? e : o;
    }

    const char* data(int j) const { return s_ + off(j); }

    size_t length(int j) const {
        size_t o = off(j), e = end(j);
        return (e >= o) ? (e - o) : 0;
    }

    // Replicate spline[jj][0] / spline[jj][2] reads on the ORIGINAL's split3() raw
    // pointers into the getline() line. Two facts drive the semantics:
    //  1. std::getline strips the trailing '\n' (hts_getline keeps it in our
    //     buffer), so the original's line ends one byte earlier; a '\r' from a
    //     CRLF file survives (the original never strips it on the VCF path).
    //  2. split3 does NOT null-terminate tokens: spline[jj] points into the
    //     line and [0]/[2] read RAW. A non-last field therefore reads past its
    //     own text into the '\t' and the next field (haploid "0" -> ch2 = next
    //     sample's first char; empty field -> ch0 = '\t'); only the last field
    //     is bounded by the line end, where reads hit '\0'. Fields of >= 3
    //     chars read their own token[0]/token[2] and are unchanged.
    size_t orig_len() const {
        size_t e = len_;
        if (e > 0 && s_[e - 1] == '\n') e--;
        return e;  // length of the line as getline() delivered it to the original
    }

    char ch0(int j) const {
        size_t o = off(j);
        return (o < orig_len()) ? s_[o] : '\0';
    }

    char ch2(int j) const {
        size_t o = off(j);
        return (o + 2 < orig_len()) ? s_[o + 2] : '\0';
    }

    std::string field(int j) const {
        size_t o = off(j), e = end(j);
        return std::string(s_ + o, (e >= o) ? (e - o) : 0);
    }

  private:
    const char* s_;
    size_t len_;
    std::vector<size_t>& offs_;
};

// Raw per-contig accumulation buffers (SoA, direct build, no per-site copies).
struct ChrBuf {
    std::vector<int64_t> pos;
    std::vector<uint64_t> minor;
    std::vector<uint64_t> miss;
};

// Empty column set used when no SubPop groups are given (all samples).
const std::vector<int> kEmptyCols;

// Builds the encoded genotype for one site into `enc` and appends the packed
// words to each group buffer. Replicates original Read_VCF_IN / Read_VCF_IN_Phase
// filters and encoding exactly. When SubPop groups are given, filters + packing
// run per group over that group's sample columns (like the original
// Read_SubPopVCF_IN*); with no groups (`subs` empty) a single all-samples
// buffer is filled. The indel filter is site-level (REF/ALT), applied once per
// line before any group. `offs` is a reusable tab-offset scratch. Returns true
// if the site was kept for buffer 0 (all-samples / first group).
// When `rset`/`region_label` are non-null (fast-path -L region tasks) the site
// must be OWNED by the given region: RegionSet::find assigns a site in
// overlapping regions to the region with the largest start, exactly like the
// fallback read_vcf path. Without this, a site shared by overlapping regions
// would be processed once per containing task and its pairs double-counted.
bool process_data_line(const char* line, size_t line_len,
                       const std::vector<SubPop>& subs, int nsamples,
                       bool phased, const Options& opt,
                       std::vector<ChrBuf*>& bufs, ReadStats& stats,
                       std::vector<size_t>& offs,
                       const RegionSet* rset = nullptr,
                       const std::string* region_label = nullptr) {
    LineTokens tk(line, line_len, offs);
    if (rset && region_label) {
        if (rset->find(tk.data(0), tk.length(0), strtoll(tk.data(1), nullptr, 10)) !=
            *region_label) {
            return false;  // belongs to a different (overlapping) region
        }
    }
    stats.sites_in++;

    // Indel filter: Base_len = max(REF.len, ALT allele lens); a site-level
    // property, independent of the sample subset.
    int base_len = static_cast<int>(tk.length(3));
    const char* alt = tk.data(4);
    const size_t alt_len = tk.length(4);
    size_t start = 0;
    for (size_t i = 0; i <= alt_len; ++i) {
        if (i == alt_len || alt[i] == ',') {
            int al = static_cast<int>(i - start);
            if (al > base_len) base_len = al;
            start = i + 1;
        }
    }
    if (base_len > 1) {
        stats.skip_indel++;
        return false;
    }

    bool kept0 = false;
    for (size_t g = 0; g < bufs.size(); ++g) {
        const std::vector<int>& sample_cols =
            subs.empty() ? kEmptyCols : subs[g].sample_cols;
        const int n_eff = subs.empty() ? nsamples : subs[g].n_samples;
        const size_t nwords = (2 * static_cast<size_t>(n_eff) + 63) / 64;
        ChrBuf& buf = *bufs[g];

        // Count alleles as raw GT chars ('0','1','.', ...) over the used
        // samples. The 256-int count array is thread_local so it is
        // zero-initialized once per thread instead of memset per site; only
        // the touched slots are reset at the top of each site/group iteration
        // (all `continue` exits below leave the array dirty, so the reset is
        // unconditional here). The later 0..255 ascending scan for the
        // best/second-best allele is unchanged, so behavior is identical.
        static thread_local std::array<int, 256> count{};
        static thread_local unsigned char touched[256];
        static thread_local int ntouched = 0;
        for (int t = 0; t < ntouched; ++t) count[touched[t]] = 0;
        ntouched = 0;
        int het_count = 0;
        int miss_count = 0;
        auto add = [&](char ch) {
            const int c = static_cast<unsigned char>(ch);
            if (count[c]++ == 0) touched[ntouched++] = static_cast<unsigned char>(c);
        };
        for (int k = 0; k < n_eff; ++k) {
            const int j = sample_cols.empty() ? (9 + k) : sample_cols[k];
            char c0 = tk.ch0(j);
            if (c0 == '.') {
                miss_count++;
            } else {
                char c2 = tk.ch2(j);
                if (c0 != c2) het_count++;
                add(c0);
                add(c2);
            }
        }

        if ((miss_count * 1.0 / n_eff) > opt.miss) continue;
        if ((het_count * 1.0 / n_eff) > opt.het) continue;

        // best / second-best allele (skip 'N'), replicate map<char,int> iteration.
        char best_base = 'N';
        char sed_base = 'N';
        int max_c = 0;
        int sed_c = 0;
        int base_count = 0;
        for (int c = 0; c < 256; ++c) {
            if (c == 'N') continue;
            int cnt = count[c];
            if (cnt == 0) continue;
            if (cnt > max_c) {
                sed_c = max_c;
                sed_base = best_base;
                max_c = cnt;
                best_base = static_cast<char>(c);
            } else if (cnt >= sed_c) {
                sed_c = cnt;
                sed_base = static_cast<char>(c);
            }
            base_count++;
        }
        if (base_count == 1 || base_count > 2) {
            stats.skip_non_biallelic++;
            continue;
        }
        if ((sed_c * 1.0 / (sed_c + max_c)) < opt.maf) continue;

        // Encode (0/1/2) per haplotype, packing straight into the per-contig
        // flat buffers (no per-site temporaries).
        const size_t wbase = buf.minor.size();
        buf.minor.resize(wbase + nwords);
        buf.miss.resize(wbase + nwords);
        for (int k = 0; k < n_eff; ++k) {
            const int j = sample_cols.empty() ? (9 + k) : sample_cols[k];
            char c0 = tk.ch0(j);
            char a, b;
            if (c0 == '.') {
                a = 'N';
                b = 'N';
            } else {
                char c2 = tk.ch2(j);
                if (phased) {
                    a = c0;
                    b = c2;
                } else {
                    if (c0 != c2) {
                        a = best_base;
                        b = sed_base;
                    } else {
                        a = c0;
                        b = c2;
                    }
                }
            }
            auto enc = [&](char g2) -> uint8_t {
                if (g2 == best_base) return 0;
                if (g2 == sed_base) return 1;
                return 2;
            };
            const size_t h = 2 * static_cast<size_t>(k);
            for (size_t kk = 0; kk < 2; ++kk) {
                uint8_t v = enc(kk == 0 ? a : b);
                size_t w = (h + kk) / 64;
                size_t bit = (h + kk) % 64;
                if (v == 1) buf.minor[wbase + w] |= (uint64_t(1) << bit);
                else if (v == 2) buf.miss[wbase + w] |= (uint64_t(1) << bit);
            }
        }

        // Position via strtoll on the kstring buffer (null-terminated; the token
        // ends at a '\t', so no explicit terminator is needed).
        int64_t pos = strtoll(tk.data(1), nullptr, 10);

        buf.pos.push_back(pos);
        stats.sites_kept++;
        if (g == 0) kept0 = true;
    }
    return kept0;
}

std::vector<ContigData> build_contigs(std::map<std::string, ChrBuf>& raw,
                                      int nsamples) {
    const size_t nwords = (2 * static_cast<size_t>(nsamples) + 63) / 64;
    const size_t nhaplo = 2 * static_cast<size_t>(nsamples);

    std::vector<ContigData> out;
    out.reserve(raw.size());
    // Each converted ChrBuf is erased from the map as soon as its arena copy is
    // done, so the flat accumulation buffers (potentially the largest allocation
    // for a big dataset) do not coexist with the built SnpData arenas: peak RSS
    // drops from ~2x to ~1x the packed size in the single-pass fallback.
    for (auto it = raw.begin(); it != raw.end();) {
        const std::string name = it->first;
        ChrBuf& buf = it->second;
        const size_t n = buf.pos.size();
        ContigData cd;
        cd.name = name;
        cd.nwords = nwords;
        cd.n_haplotypes = nhaplo;
        cd.snps.reserve(n);

        // Stable sort by pos (replicates std::map<pos,...> order); keep the
        // first occurrence of each position (replicates map insert semantics).
        // VCF is normally already sorted per contig, so skip the sort in that
        // common case (identity ord keeps the byte-identical result).
        std::vector<size_t> ord(n);
        for (size_t i = 0; i < n; ++i) ord[i] = i;
        bool sorted = true;
        for (size_t i = 1; i < n; ++i) {
            if (buf.pos[i - 1] > buf.pos[i]) {
                sorted = false;
                break;
            }
        }
        if (!sorted) {
            std::stable_sort(ord.begin(), ord.end(),
                             [&](size_t x, size_t y) { return buf.pos[x] < buf.pos[y]; });
        }

        for (size_t k = 0; k < n; ++k) {
            const size_t i = ord[k];
            if (k > 0 && buf.pos[ord[k - 1]] == buf.pos[i]) continue;  // dedupe
            SnpData sd;
            sd.pos = buf.pos[i];
            sd.minor.assign(cd.arena, &buf.minor[i * nwords], nwords);
            sd.miss.assign(cd.arena, &buf.miss[i * nwords], nwords);
            cd.snps.push_back(std::move(sd));
        }
        out.push_back(std::move(cd));
        it = raw.erase(it);
    }
    return out;
}

}  // namespace

int read_vcf(const std::string& path, const Options& opt,
             std::vector<std::vector<ContigData>>& outs, ReadStats& stats,
             const std::vector<SubPop>& subs) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "open VCF File IN File error: " << path << std::endl;
        return 0;
    }

    // Parallel BGZF decompression for compressed inputs (transparent: records
    // come out in identical order, so byte-compatibility is unaffected).
    if (opt.threads > 1 && hts_get_bgzfp(fp)) {
        hts_set_threads(fp, opt.threads);
    }

    // One buffer set per SubPop group (or a single all-samples set).
    const size_t nbufs = subs.empty() ? 1 : subs.size();
    std::vector<std::map<std::string, ChrBuf>> raws(nbufs);
    int nsamples = -1;
    bool phased = false;
    bool saw_data = false;
    RegionSet rset(opt.regions);
    std::vector<size_t> offs;  // reusable tab-offset scratch (one per reader)
    // Reused across data lines (targets re-pointed per line; capacity kept).
    std::vector<ChrBuf*> gbufs;
    gbufs.reserve(nbufs);

    // OutFilterSNP output: the passing data lines in original file order,
    // written as a plain gzip (matches the original ogzstream output). Only in
    // the full-sample readers (SubPop+OutFilterSNP writes nothing, as the
    // original). The main() forces the fallback path when this is requested.
    gzFile filt = nullptr;
    std::string filt_path;
    if (opt.out_filter_snp && subs.empty()) {
        filt_path = strip_stat_prefix(opt.out_stat) + ".vcf.filter.gz";
        filt = gzopen(filt_path.c_str(), "wb");
        if (!filt) {
            std::cerr << "open Filter SNP OUT File error: " << filt_path << std::endl;
            hts_close(fp);
            return 0;
        }
    }

    kstring_t line = {0, 0, 0};
    int ret;
    while ((ret = hts_getline(fp, '\n', &line)) >= 0) {
        if (line.l == 0) continue;
        if (line.s[0] == '#') {
            if (nsamples >= 0) continue;  // header already consumed
            if (line.l > 1 && line.s[1] == '#') continue;
            std::string l(line.s, line.l);
            std::vector<std::string> toks = split_ws(l);
            if (toks.empty() || toks[0] != "#CHROM") continue;
            nsamples = static_cast<int>(toks.size()) - 9;
            if (nsamples < 0) nsamples = 0;
            continue;
        }
        if (nsamples < 0) {
            std::cerr << "VCF Header something wrong, can't find sample info before site info" << std::endl;
            hts_close(fp);
            return 0;
        }
        if (!saw_data) {
            LineTokens tk(line.s, line.l, offs);
            size_t o = tk.off(9);
            phased = (o + 1 < line.l) && line.s[o + 1] == '|';
            saw_data = true;
        }
        LineTokens tk(line.s, line.l, offs);
        if (!opt.regions.empty()) {
            // Region restrict: keep only sites inside some -L region. Each
            // region is an independent accumulation unit, so cross-region
            // pairs are excluded. Sites outside all regions are dropped here
            // (before the MAF/het/miss filters), like bcftools -r. The chrom
            // name is matched length-bounded against the RegionSet (no
            // temporary string); empty/blank records drop out like the map key
            // case below.
            int64_t pos = strtoll(tk.data(1), nullptr, 10);
            const std::string& label = rset.find(tk.data(0), tk.length(0), pos);
            if (label.empty()) continue;
            gbufs.clear();
            for (size_t g = 0; g < nbufs; ++g) gbufs.push_back(&raws[g][label]);
            if (process_data_line(line.s, line.l, subs, nsamples, phased, opt,
                                  gbufs, stats, offs)) {
                if (filt) gzwrite(filt, line.s, static_cast<unsigned>(line.l));
                if (filt) gzputc(filt, '\n');
            }
        } else {
            std::string chrom = tk.field(0);
            if (chrom.empty()) continue;
            gbufs.clear();
            for (size_t g = 0; g < nbufs; ++g) gbufs.push_back(&raws[g][chrom]);
            if (process_data_line(line.s, line.l, subs, nsamples, phased, opt,
                                  gbufs, stats, offs)) {
                if (filt) gzwrite(filt, line.s, static_cast<unsigned>(line.l));
                if (filt) gzputc(filt, '\n');
            }
        }
    }
    if (filt) gzclose(filt);
    hts_close(fp);

    outs.clear();
    outs.reserve(nbufs);
    for (size_t g = 0; g < nbufs; ++g) {
        const int n_eff = subs.empty() ? nsamples : subs[g].n_samples;
        outs.push_back(build_contigs(raws[g], n_eff));
    }
    return 1;
}

bool compute_bins_indexed(const std::string& path, const Options& opt,
                          int max_dist_bp, std::vector<AccumSpec>& outs,
                          ReadStats& stats, const std::vector<SubPop>& subs) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const bool timing = getenv("PLD2_TIMING") != nullptr;
    const auto ms = [](const Clock::time_point& t) {
        return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
    };

    const size_t nbufs = subs.empty() ? 1 : subs.size();
    const int out_type = outs[0].out_type;
    const size_t nbins = static_cast<size_t>(max_dist_bp) + 1;
    const bool hist_mode = (out_type == 4 || out_type == 5);
    for (AccumSpec& out : outs) {
        out.bins.assign(nbins, Bin{});
        if (hist_mode) {
            out.hist_rr.assign(nbins * 101, 0);
            out.hist_d.assign((out_type == 4) ? nbins * 101 : 0, 0);
        }
        if (out_type == 3 || out_type == 6 || out_type == 7 || out_type == 8) {
            out.ld = std::make_unique<LdWriter>(out_type);
        }
    }

    // ---- prelude: header -> nsamples, contig list, phase (replicates old tool) ----
    htsFile* pre = hts_open(path.c_str(), "r");
    if (!pre) {
        std::cerr << "open VCF File IN File error: " << path << std::endl;
        return false;
    }
    int nsamples = -1;
    bool phased = false;
    bool got_phase = false;
    std::vector<std::string> contigs;
    std::vector<int64_t> contig_lens;  // -1 when no length= in ##contig
    kstring_t line = {0, 0, 0};
    while (hts_getline(pre, '\n', &line) >= 0) {
        if (line.l == 0) continue;
        std::string l(line.s, line.l);
        if (l[0] == '#') {
            if (l.size() > 1 && l[1] == '#') {
                const std::string tag = "##contig=<ID=";
                size_t p = l.find(tag);
                if (p != std::string::npos) {
                    size_t start = p + tag.size();
                    size_t q = l.find_first_of(",>", start);
                    if (q != std::string::npos) {
                        contigs.push_back(l.substr(start, q - start));
                        int64_t clen = -1;
                        size_t lp = l.find("length=", q);
                        if (lp != std::string::npos) {
                            size_t le = lp + 7;
                            while (le < l.size() && isdigit(
                                                       static_cast<unsigned char>(l[le])))
                                ++le;
                            if (le > lp + 7) {
                                clen = std::strtoll(l.c_str() + lp + 7, nullptr, 10);
                            }
                        }
                        contig_lens.push_back(clen);
                    }
                }
                continue;
            }
            std::vector<std::string> toks = split_ws(l);
            if (toks.empty() || toks[0] != "#CHROM") continue;
            nsamples = static_cast<int>(toks.size()) - 9;
            if (nsamples < 0) nsamples = 0;
            continue;
        }
        if (!got_phase) {
            std::vector<size_t> offs;
            LineTokens tk(line.s, line.l, offs);
            size_t o = tk.off(9);
            phased = (o + 1 < line.l) && line.s[o + 1] == '|';
            got_phase = true;
        }
        break;  // only the first data line is needed for phase detection
    }
    hts_close(pre);

    if (nsamples < 0 || !got_phase || contigs.empty()) {
        if (timing) std::cerr << "[timing] fastpath declined (header/index conditions)\n";
        return false;
    }

    tbx_t* tbx = tbx_index_load(path.c_str());
    if (!tbx) {
        if (timing) std::cerr << "[timing] fastpath declined (no .tbi/.csi)\n";
        return false;
    }

    // Per-group sample counts / word counts / haplotype counts.
    std::vector<size_t> n_eff(nbufs), nwords(nbufs), nhaplo(nbufs);
    for (size_t g = 0; g < nbufs; ++g) {
        n_eff[g] = subs.empty() ? static_cast<size_t>(nsamples)
                                : static_cast<size_t>(subs[g].n_samples);
        nwords[g] = (2 * n_eff[g] + 63) / 64;
        nhaplo[g] = 2 * n_eff[g];
        outs[g].nhaplo = nhaplo[g];
    }

    // Task contigs come from the INDEX (tbx_seqnames), not the header: a VCF may
    // carry data on contigs that have no ##contig=<ID=...> line, and a
    // header-only task list would silently drop those contigs' pairs (the
    // single-pass fallback reads every chrom). The ##contig lines still provide
    // the declared lengths that drive position sharding; index-only contigs get
    // no length (single task each). Contig names that would break tabix
    // region-query parsing (':' splits chrom/pos) make tbx_itr_querys return
    // nothing for the whole contig, so decline the fast path for those (the
    // map-keyed fallback handles such names fine).
    std::vector<std::string> task_chroms;
    std::vector<int64_t> task_lens;
    {
        std::map<std::string, int64_t> len_by_name;
        for (size_t i = 0; i < contigs.size(); ++i) len_by_name[contigs[i]] = contig_lens[i];
        const auto unsafe = [](const std::string& n) {
            return n.empty() || n.find(':') != std::string::npos;
        };
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (unsafe(contigs[i])) {
                if (timing) std::cerr << "[timing] fastpath declined (contig name unsafe for tabix query)\n";
                return false;
            }
        }
        int nidx = 0;
        const char** seqs = tbx_seqnames(tbx, &nidx);
        for (int i = 0; i < nidx; ++i) {
            const std::string n(seqs[i]);
            if (unsafe(n)) {
                std::free(seqs);
                if (timing) std::cerr << "[timing] fastpath declined (contig name unsafe for tabix query)\n";
                return false;
            }
        }
        // Header-declared contigs that actually have data (header order, using
        // their declared lengths), then any index-only contigs (no declared
        // length -> single task each).
        std::set<std::string> in_index(seqs, seqs + nidx);
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (in_index.count(contigs[i])) {
                task_chroms.push_back(contigs[i]);
                task_lens.push_back(contig_lens[i]);
            }
        }
        for (int i = 0; i < nidx; ++i) {
            const std::string n(seqs[i]);
            if (!len_by_name.count(n)) {
                task_chroms.push_back(n);
                task_lens.push_back(-1);
            }
        }
        std::free(seqs);
    }

    // Task units: with -L each region is an independent accumulation unit
    // (cross-region pairs excluded); otherwise one task per contig, or -- when
    // more threads than contigs -- position shards of the larger contigs so the
    // worker pool stays saturated on low-chromosome datasets.
    struct FastTask {
        std::string query;       // tabix region string
        std::string row_label;   // cd.name written into emitted rows
        std::string part_label;  // LdWriter part label (sort key for .LD.gz)
        int64_t core_hi = -1;    // exclusive left-site bound; -1 = whole buffer
    };
    // Upper position bound for the open-ended last shard (see below): the VCF's
    // declared ##contig length is only an upper bound we can trust for query
    // windows when it is not the *last* shard -- the final shard of each contig
    // must reach the end of the chromosome so sites beyond the declared length
    // (a mis-annotated / inconsistent reference) are still fetched and owned.
    constexpr int64_t kMaxPos = INT64_MAX - 1;
    std::vector<FastTask> tasks;
    if (!opt.regions.empty()) {
        for (const Region& r : opt.regions) {
            tasks.push_back(FastTask{region_query_string(r), r.label, r.label, -1});
        }
    } else {
        const size_t ncontigs = task_chroms.size();
        bool any_split = false;
        std::vector<size_t> kshards(ncontigs, 1);
        if (static_cast<size_t>(opt.threads) > ncontigs) {
            // Target span so sharding fills the thread budget; length-weighted.
            int64_t total_len = 0;
            for (int64_t clen : task_lens)
                if (clen > 0) total_len += clen;
            if (total_len > 0) {
                const int64_t span = total_len / opt.threads + 1;
                const int64_t min_span = 1000000;  // never split below ~1 Mb
                for (size_t i = 0; i < ncontigs; ++i) {
                    const int64_t clen = task_lens[i];
                    if (clen <= 0) continue;
                    int64_t k = (clen + span / 2) / span;
                    if (k < 1) k = 1;
                    if (clen / min_span < k) k = clen / min_span;
                    if (k < 1) k = 1;
                    kshards[i] = static_cast<size_t>(k);
                }
            }
            for (size_t k : kshards)
                if (k > 1) any_split = true;
        }
        if (any_split) {
            for (size_t i = 0; i < ncontigs; ++i) {
                const std::string& chr = task_chroms[i];
                const size_t k = kshards[i];
                const int64_t clen = task_lens[i];
                if (k == 1 || clen <= 0) {
                    tasks.push_back(FastTask{chr, chr, chr, -1});
                    continue;
                }
                const int64_t span_i = (clen + static_cast<int64_t>(k) - 1) /
                                       static_cast<int64_t>(k);
                char pad[32];
                for (size_t s = 0; s < k; ++s) {
                    const bool last = (s + 1 == k);
                    const int64_t lo = 1 + static_cast<int64_t>(s) * span_i;
                    // The last shard owns every site from `lo` to the end of the
                    // chromosome: bound it by kMaxPos (not clen+1) so sites past
                    // the declared ##contig length are fetched and left-partners
                    // (their pairs must be counted exactly once, like the
                    // unsharded fallback). Earlier shards stay bounded by clen+1.
                    const int64_t hi_excl =
                        last ? kMaxPos
                             : std::min<int64_t>(clen + 1,
                                                 1 + static_cast<int64_t>(s + 1) * span_i);
                    std::snprintf(pad, sizeof(pad), "%06zu", s);
                    const std::string plabel = chr + '\x01' + pad;
                    const std::string query =
                        last ? (chr + ":" + std::to_string(lo) + "-" +
                                std::to_string(kMaxPos))
                             : (chr + ":" + std::to_string(lo) + "-" +
                                std::to_string(hi_excl + max_dist_bp));
                    tasks.push_back(FastTask{std::move(query), chr, plabel, hi_excl});
                }
            }
        } else {
            for (size_t i = 0; i < ncontigs; ++i) {
                tasks.push_back(FastTask{task_chroms[i], task_chroms[i], task_chroms[i], -1});
            }
        }
    }
    if (tasks.empty()) {
        if (timing) std::cerr << "[timing] fastpath declined (no tasks)\n";
        return false;
    }
    size_t nworkers = std::min<size_t>(static_cast<size_t>(opt.threads), tasks.size());
    if (hist_mode && nworkers > 4) {
        // OutType 4/5 keeps a (MaxDist+1)*101 int histogram per worker
        // (~121 MB/worker at MaxDist=300 kb). Cap the worker count to bound
        // memory; integer sums merge in any order, so output is byte-identical.
        std::cerr << "OutType 4/5 histogram mode: capping accumulation workers "
                     "to 4 (histogram is ~121 MB/worker at MaxDist=300 kb)\n";
        nworkers = 4;
    }
    if (timing) std::cerr << "[timing] fastpath prelude=" << ms(t0) << "ms  nworkers=" << nworkers
                          << " nsamples=" << nsamples << " phased=" << phased
                          << " tasks=" << tasks.size() << "\n";

    // Per-worker partials, one AccumSpec per group.
    std::vector<std::vector<AccumSpec>> partial(nworkers);
    for (auto& pw : partial) pw.resize(nbufs);
    for (auto& pw : partial) {
        for (size_t g = 0; g < nbufs; ++g) {
            AccumSpec& p = pw[g];
            p.out_type = out_type;
            p.nhaplo = nhaplo[g];
            p.bins.assign(nbins, Bin{});
            if (hist_mode) {
                p.hist_rr.assign(nbins * 101, 0);
                if (out_type == 4) p.hist_d.assign(nbins * 101, 0);
            }
            if (out_type == 3 || out_type == 6 || out_type == 7 || out_type == 8) {
                p.ld = std::make_unique<LdWriter>(out_type);
            }
        }
    }
    std::vector<ReadStats> wstats(nworkers);
    std::atomic<size_t> next{0};
    std::atomic<bool> worker_fail{false};

    // -L region ownership (fastpath): each region is a separate tabix query, so
    // without a per-site ownership check, sites inside overlapping regions
    // would be processed once per containing task (pairs double-counted vs the
    // fallback). Build the fallback's RegionSet once and pass it to each task.
    RegionSet rset(opt.regions);
    const bool region_filter = !opt.regions.empty();

    const auto worker = [&](size_t wid) {
        htsFile* fp = hts_open(path.c_str(), "r");
        if (!fp) {
            worker_fail.store(true);
            return;
        }
        // Reusable per-task accumulation buffers: kept across tasks so shard
        // processing does not reallocate ChrBuf storage for every task.
        std::vector<ChrBuf> bufs(nbufs);
        // Group buffer pointers are stable across line/task boundaries (clear()
        // never invalidates), so fill once instead of per data line.
        std::vector<ChrBuf*> gbufs;
        gbufs.reserve(nbufs);
        for (size_t g = 0; g < nbufs; ++g) gbufs.push_back(&bufs[g]);
        for (;;) {
            const size_t ci = next.fetch_add(1);
            if (ci >= tasks.size()) break;
            hts_itr_t* itr = tbx_itr_querys(tbx, tasks[ci].query.c_str());
            if (!itr) {
                // A failed query would silently drop the whole task's pairs;
                // treat it as a worker failure so main() falls back to the
                // single-pass reader instead of writing incomplete output.
                worker_fail.store(true);
                break;
            }
            kstring_t ln = {0, 0, 0};
            std::vector<size_t> offs;  // reusable tab-offset scratch (per worker)
            int ret;
            while ((ret = tbx_itr_next(fp, tbx, itr, &ln)) >= 0) {
                if (ln.l == 0) continue;
                process_data_line(ln.s, ln.l, subs, nsamples, phased, opt, gbufs,
                                  wstats[wid], offs,
                                  region_filter ? &rset : nullptr,
                                  region_filter ? &tasks[ci].row_label : nullptr);
            }
            hts_itr_destroy(itr);

            // Indexed files are sorted per contig: records arrive in ascending
            // pos order, so dedupe is adjacent (keep first, like the map path).
            for (size_t g = 0; g < nbufs; ++g) {
                ContigData cd;
                cd.name = tasks[ci].row_label;
                cd.nwords = nwords[g];
                cd.n_haplotypes = nhaplo[g];
                const size_t n = bufs[g].pos.size();
                cd.snps.reserve(n);
                int64_t last = -1;
                bool have_last = false;
                for (size_t i = 0; i < n; ++i) {
                    if (have_last && bufs[g].pos[i] == last) continue;
                    last = bufs[g].pos[i];
                    have_last = true;
                    SnpData sd;
                    sd.pos = bufs[g].pos[i];
                    sd.minor.assign(cd.arena, &bufs[g].minor[i * nwords[g]], nwords[g]);
                    sd.miss.assign(cd.arena, &bufs[g].miss[i * nwords[g]], nwords[g]);
                    cd.snps.push_back(std::move(sd));
                }
                // Shard left-site core: pairs are owned by the shard whose core
                // window holds the left site (deduped, positions ascending).
                size_t core_n = cd.snps.size();
                if (tasks[ci].core_hi >= 0) {
                    size_t lo = 0, hi = cd.snps.size();
                    while (lo < hi) {
                        const size_t mid = lo + (hi - lo) / 2;
                        if (cd.snps[mid].pos < tasks[ci].core_hi) lo = mid + 1;
                        else hi = mid;
                    }
                    core_n = lo;
                }
                if (core_n == 0) continue;
                accumulate_contig(cd, max_dist_bp, partial[wid][g], core_n,
                                  tasks[ci].part_label);
            }
            // Reuse the accumulation buffers for the next task (capacity kept).
            for (ChrBuf& b : bufs) {
                b.pos.clear();
                b.minor.clear();
                b.miss.clear();
            }
        }
        hts_close(fp);
    };

    std::vector<std::thread> pool;
    pool.reserve(nworkers - 1);
    for (size_t w = 1; w < nworkers; ++w) pool.emplace_back(worker, w);
    worker(0);
    for (auto& t : pool) t.join();
    tbx_destroy(tbx);

    if (worker_fail.load()) {
        if (timing) std::cerr << "[timing] fastpath failed (worker open)\n";
        return false;
    }

    // Merge partials in worker order per group (deterministic; double sums
    // ordered per worker, see ldkernel.h note).
    for (size_t g = 0; g < nbufs; ++g) {
        AccumSpec& out = outs[g];
        for (size_t w = 0; w < nworkers; ++w) {
            for (size_t d = 0; d < nbins; ++d) {
                out.bins[d].count += partial[w][g].bins[d].count;
                out.bins[d].sumRR += partial[w][g].bins[d].sumRR;
                out.bins[d].sumD += partial[w][g].bins[d].sumD;
            }
            if (hist_mode) {
                for (size_t k = 0; k < nbins * 101; ++k) {
                    out.hist_rr[k] += partial[w][g].hist_rr[k];
                    if (out_type == 4) out.hist_d[k] += partial[w][g].hist_d[k];
                }
            }
            if (out.ld && partial[w][g].ld) out.ld->absorb(*partial[w][g].ld);
        }
    }
    for (size_t w = 0; w < nworkers; ++w) {
        stats.sites_in += wstats[w].sites_in;
        stats.sites_kept += wstats[w].sites_kept;
        stats.skip_indel += wstats[w].skip_indel;
        stats.skip_non_biallelic += wstats[w].skip_non_biallelic;
    }
    if (timing) std::cerr << "[timing] fastpath total=" << ms(t0) << "ms\n";
    return true;
}

}  // namespace pld2