#include "ehh.h"

#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace pld2 {

namespace {

// Replicates the original split(): splits on any run of the delimiter,
// skipping empty fields.
void split_colon(const std::string& s, std::vector<std::string>& out) {
    const std::string delim = ":";
    std::string::size_type last = s.find_first_not_of(delim, 0);
    std::string::size_type pos = s.find_first_of(delim, last);
    while (pos != std::string::npos || last != std::string::npos) {
        out.push_back(s.substr(last, pos - last));
        last = s.find_first_not_of(delim, pos);
        pos = s.find_first_of(delim, last);
    }
}

// Single-quotes a string for use in a system() shell command (embedded single
// quotes become the POSIX '\'' splice). Handles spaces and shell metacharacters
// in user-supplied output prefixes.
std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out.push_back(c);
        }
    }
    out += "'";
    return out;
}

// Escapes a string for embedding inside a double-quoted R string literal
// (backslash and double-quote only).
std::string r_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == '\\' || c == '"') out.push_back('\\');
        out.push_back(c);
    }
    return out;
}

struct EhhSite {
    int64_t pos;
    std::string haps;  // Asize chars: '0'/'1' (no-missing sites only)
};

// Groups the current haplotype strings, replicating the original's
// std::map<string,int> count + fengzi/fengzi_A accumulation (integer sums, so
// the container/order do not matter).
void count_groups(const std::vector<std::string>& haplotypes, int64_t& fengzi,
                  int64_t& fengzi_A) {
    std::unordered_map<std::string, int> cnt;
    for (const std::string& h : haplotypes) cnt[h]++;
    fengzi = 0;
    fengzi_A = 0;
    for (const auto& kv : cnt) {
        int64_t c = static_cast<int64_t>(kv.second) * (kv.second - 1);
        fengzi += c;
        if (kv.first[0] == '0') fengzi_A += c;
    }
}

}  // namespace

bool parse_ehh(const std::string& spec, EhhParams& out) {
    std::vector<std::string> tok;
    split_colon(spec, tok);
    if (tok.size() < 2) return false;
    out.chr = tok[0];
    out.start_site = std::atoi(tok[1].c_str());
    return true;
}

int run_ehh(const std::vector<ContigData>& contigs, const EhhParams& p,
            int max_dist_bp, int tf, const std::string& prefix) {
    const std::string out_path = prefix + ".ehh.gz";
    gzFile OUT = gzopen(out_path.c_str(), "wb");
    if (!OUT) {
        std::cerr << "open OUT File error: " << out_path << std::endl;
        return 0;
    }

    const ContigData* cd = nullptr;
    for (const ContigData& c : contigs) {
        if (c.name == p.chr) {
            cd = &c;
            break;
        }
    }
    if (!cd || cd->snps.empty()) {
        std::cerr << "\t\tInPut Para -EHH  chromosome [" << p.chr
                  << "]  can't be found in the SNP dataset\n";
        gzclose(OUT);
        return 0;
    }

    const int64_t begin = p.start_site - max_dist_bp < 0 ? 0 : p.start_site - max_dist_bp;
    const int64_t end = p.start_site + max_dist_bp;
    const size_t n_haps = cd->n_haplotypes;
    const int64_t Asize = static_cast<int64_t>(n_haps);

    std::vector<EhhSite> sites;
    int TotalSNP = 0;
    for (const SnpData& snp : cd->snps) {
        if (snp.pos < begin) continue;
        if (snp.pos > end) break;
        TotalSNP++;
        bool no_missing = true;
        for (size_t h = 0; h < n_haps; h++) {
            if (snp.miss[h >> 6] >> (h & 63) & 1ULL) {
                no_missing = false;
                break;
            }
        }
        if (!no_missing) continue;
        std::string haps(n_haps, '0');
        for (size_t h = 0; h < n_haps; h++) {
            haps[h] = static_cast<char>('0' + (snp.minor[h >> 6] >> (h & 63) & 1ULL));
        }
        sites.push_back({snp.pos, std::move(haps)});
    }

    std::cout << "##Start EHH region :" << p.chr << " " << begin << " " << end
              << "; In This Region TotalSNP Number is " << TotalSNP
              << ",No Missing SNP Site Number is " << sites.size() << std::endl;

    const int64_t count = static_cast<int64_t>(sites.size());
    if (count < 8) {
        std::cerr << "EHH should be No Missing SNP Number  > 8 \n";
        gzclose(OUT);
        return 0;
    }
    if (count > 168888) {
        std::cerr << "Warning: EHH Region SNP Number too much,you may use the small region or more stringent conditions to filter the SNP\n";
    }

    const int64_t fengmu = Asize * (Asize - 1);

    size_t core = sites.size();
    for (size_t i = 0; i < sites.size(); i++) {
        if (sites[i].pos <= p.start_site) {
            core = i;
        } else {
            break;
        }
    }
    if (core == 0) core = 1;
    if (core == sites.size()) core -= 2;

    std::unordered_map<int64_t, std::pair<double, double>> result;  // pos -> (D, RR)

    std::vector<std::string> haplotypes(n_haps);
    int64_t CalTmp = 0;
    for (size_t h = 0; h < n_haps; h++) {
        haplotypes[h].push_back(sites[core].haps[h]);
        if (sites[core].haps[h] == '0') CalTmp++;
    }
    double ehh_rr = (static_cast<double>(CalTmp * (CalTmp - 1))) / fengmu;
    CalTmp = Asize - CalTmp;
    double ehh_d = (static_cast<double>(CalTmp * (CalTmp - 1))) / fengmu + ehh_rr;
    result[sites[core].pos] = {ehh_d, ehh_rr};

    std::cout << "##Begin Cal EHH...\n";

    // Left loop: extend core haplotypes leftwards, appending one char per site.
    size_t key2 = core - 1;
    for (;; key2--) {
        for (size_t h = 0; h < n_haps; h++) haplotypes[h].push_back(sites[key2].haps[h]);
        if (haplotypes[0].empty() || haplotypes[0].size() < 2) continue;
        int64_t fengzi = 0, fengzi_A = 0;
        count_groups(haplotypes, fengzi, fengzi_A);
        ehh_d = static_cast<double>(fengzi) / fengmu;
        ehh_rr = static_cast<double>(fengzi_A) / fengmu;
        result[sites[key2].pos] = {ehh_d, ehh_rr};
        if (ehh_d < 0.088) break;
        if (key2 == 0) break;
    }

    // Right loop: fresh from the core, extending rightwards.
    for (size_t h = 0; h < n_haps; h++) haplotypes[h].assign(1, sites[core].haps[h]);
    for (key2 = core + 1; key2 < sites.size(); key2++) {
        for (size_t h = 0; h < n_haps; h++) haplotypes[h].push_back(sites[key2].haps[h]);
        if (haplotypes[0].empty() || haplotypes[0].size() < 2) continue;
        int64_t fengzi = 0, fengzi_A = 0;
        count_groups(haplotypes, fengzi, fengzi_A);
        ehh_d = static_cast<double>(fengzi) / fengmu;
        ehh_rr = static_cast<double>(fengzi_A) / fengmu;
        result[sites[key2].pos] = {ehh_d, ehh_rr};
        if (ehh_d < 0.088) break;
    }

    std::string out;
    out.reserve(4096);
    out += "#Chr\tSite\tDist\tEHH_all\tEHH_0\tEHH_1\n";
    out += "#SNP_Number\t" + std::to_string(count) + "\n";
    char row[256];
    for (const EhhSite& s : sites) {
        const int64_t dis = s.pos - p.start_site;
        auto it = result.find(s.pos);
        if (it == result.end()) {
            out += p.chr + "\t" + std::to_string(s.pos) + "\t" +
                   std::to_string(dis) + "\t0.0000\t0.0000\t0.0000\n";
        } else {
            std::snprintf(row, sizeof(row), "%s\t%lld\t%lld\t%.4f\t%.4f\t%.4f\n",
                          p.chr.c_str(), static_cast<long long>(s.pos),
                          static_cast<long long>(dis), it->second.first,
                          it->second.second, it->second.first - it->second.second);
            out += row;
        }
    }
    gzwrite(OUT, out.data(), static_cast<unsigned>(out.size()));
    gzclose(OUT);

    // ---------- plot EHH (replicates the original's Rscript section) ----------
    // popen can fail (fork/exec failure); bail out of the R plot step then,
    // instead of fread()ing a NULL stream (the .ehh table itself is already
    // written above).
    char buf[2048] = {'\0'};
    FILE* stream = popen("which  Rscript  2> /dev/null ", "r");
    if (stream) {
        std::fread(buf, sizeof(char), sizeof(buf), stream);
        buf[sizeof(buf) - 1] = '\0';  // clamp in case `which` filled the buffer
        pclose(stream);
    }
    std::string binPath = buf;
    if (binPath.length() > 0) binPath = binPath.substr(0, binPath.length() - 1);

    const std::string OutPlotr = prefix + ".tmp.r";
    {
        FILE* f = std::fopen(OutPlotr.c_str(), "w");
        if (f) {
            const std::string e_prefix = r_escape(prefix);
            const std::string e_chr = r_escape(p.chr);
            const std::string e_out = r_escape(out_path);
            std::fprintf(f,
                "\n"
                "read.table(\"%s\")->r;\n"
                "pdf(\"%s.ehh.pdf\");\n"
                "plot(r[,3]/1000,r[,4],col=\"blue\",type=\"l\",ylab=\"EHH\",main=\"EHH decay\",,bty=\"n\",xlab=\"Distance from core region %s:%lld (Kb)\")\n"
                "dev.off();\n"
                "png(\"%s.ehh.png\");\n"
                "plot(r[,3]/1000,r[,4],col=\"blue\",type=\"l\",ylab=\"EHH\",main=\"EHH decay\",,bty=\"n\",xlab=\"Distance from core region %s:%lld (Kb)\")\n"
                "dev.off();\n"
                "\n",
                e_out.c_str(), e_prefix.c_str(), e_chr.c_str(),
                static_cast<long long>(p.start_site), e_prefix.c_str(), e_chr.c_str(),
                static_cast<long long>(p.start_site));
            std::fprintf(f, "\n");
            std::fclose(f);
        }
    }

    if (binPath.empty()) {
        std::cout << "\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"
                  << std::endl;
        std::cout << "\t\tRscript " << OutPlotr << std::endl;
    } else {
        std::string cc = shell_quote(binPath) + " " + shell_quote(OutPlotr);
        if (tf) {
            cc += " ; rm -rf " + shell_quote(OutPlotr);
        }
        std::system(cc.c_str());
    }

    std::system(("rm -rf " + shell_quote(prefix + ".stat.gz")).c_str());

    return 1;
}

}  // namespace pld2
