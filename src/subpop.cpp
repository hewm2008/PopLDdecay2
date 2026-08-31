#include "subpop.h"

#include <fstream>
#include <iostream>
#include <map>

#include <htslib/hts.h>
#include <htslib/kstring.h>

namespace pld2 {

namespace {

// Tokenization on space/tab (replicates split(line,inf," \t")).
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

// Basename with the last extension stripped ("/a/b/popA.list" -> "popA").
std::string base_name(const std::string& p) {
    const size_t slash = p.find_last_of("/\\");
    std::string b = (slash == std::string::npos) ? p : p.substr(slash + 1);
    const size_t dot = b.find_last_of('.');
    if (dot != std::string::npos && dot > 0) b = b.substr(0, dot);
    return b;
}

// Filename-safe group label (used as an output suffix).
std::string sanitize_label(const std::string& s) {
    std::string r = s;
    for (char& c : r) {
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                        (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
        if (!ok) c = '_';
    }
    if (r.empty()) r = "group";
    return r;
}

}  // namespace

bool parse_subpop_lists(const std::vector<std::string>& paths,
                        std::vector<SubGroupSpec>& out) {
    out.clear();
    if (paths.empty()) return true;
    const bool legacy = (paths.size() == 1);

    for (const std::string& p : paths) {
        std::ifstream in(p.c_str());
        if (!in) {
            std::cerr << "open Sub Group IN File error: " << p << std::endl;
            return false;
        }
        std::vector<std::vector<std::string>> rows;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
            if (line.empty()) continue;
            rows.push_back(split_ws(line));
        }

        // Column-mode detection: a file is two-column (`sample<TAB>groupid`)
        // only when EVERY data line has exactly 2 tokens. Any other layout
        // (single names per line, or the legacy multi-name-per-line style) is
        // treated as the original single-column format where every token is a
        // sample name.
        bool all_two = !rows.empty();
        for (const auto& r : rows) {
            if (r.size() != 2) {
                all_two = false;
                break;
            }
        }
        if (!all_two) {
            // Single-column: one group, label = stripped basename (or "" in
            // legacy single-file mode), names = all tokens.
            SubGroupSpec spec;
            spec.label = legacy ? "" : sanitize_label(base_name(p));
            for (const auto& r : rows) spec.names.insert(r.begin(), r.end());
            out.push_back(std::move(spec));
        } else {
            // Two-column: sample<TAB>groupid -> one group per distinct id.
            std::vector<std::string> order;
            std::map<std::string, std::set<std::string>> by_id;
            for (const auto& r : rows) {
                const std::string& id = r[1];
                if (by_id.find(id) == by_id.end()) order.push_back(id);
                by_id[id].insert(r[0]);
            }
            for (const std::string& id : order) {
                SubGroupSpec spec;
                spec.label = sanitize_label(id);
                spec.names = by_id[id];
                out.push_back(std::move(spec));
            }
        }
    }
    return true;
}

bool resolve_subpop_header(const std::vector<std::string>& header,
                           const std::vector<SubGroupSpec>& specs,
                           const char* found_msg, const char* header_kind,
                           int sample_start, std::vector<SubPop>& out) {
    out.clear();
    for (const SubGroupSpec& spec : specs) {
        // Match in header order; count occurrences per name for warnings.
        std::map<std::string, int> occ;
        for (const std::string& n : spec.names) occ[n] = 0;
        SubPop sp;
        sp.active = true;
        sp.label = spec.label;
        for (size_t j = static_cast<size_t>(sample_start); j < header.size(); ++j) {
            auto it = occ.find(header[j]);
            if (it != occ.end()) {
                sp.sample_cols.push_back(static_cast<int>(j));
                it->second++;
            }
        }
        const int n_sub = static_cast<int>(sp.sample_cols.size());
        std::cout << "the Number of subPop samples" << found_msg << " is " << n_sub << std::endl;
        if (n_sub < 3) {
            std::cerr << "sub Group Population szie is too small, SubGroup sample size: "
                      << n_sub << std::endl;
            return false;
        }
        for (const auto& kv : occ) {
            if (kv.second == 0) {
                std::cerr << "warning : Sample [ " << kv.first
                          << " ] can't be found in " << header_kind << "\n";
            } else if (kv.second > 1) {
                std::cerr << "warning : Sample [ " << kv.first
                          << " ] can be found [Repeat] in " << header_kind << "\n";
            }
        }
        sp.n_samples = n_sub;
        out.push_back(std::move(sp));
    }
    return true;
}

bool resolve_subpops(const std::string& vcf_path,
                     const std::vector<SubGroupSpec>& specs,
                     std::vector<SubPop>& out) {
    if (specs.empty()) {
        out.clear();
        return true;
    }

    // ---- VCF header -> sample columns (read once, shared by all groups) ----
    htsFile* fp = hts_open(vcf_path.c_str(), "r");
    if (!fp) {
        std::cerr << "open VCF File IN File error: " << vcf_path << std::endl;
        return false;
    }
    std::vector<std::string> header;
    kstring_t ks = {0, 0, 0};
    bool found = false;
    while (hts_getline(fp, '\n', &ks) >= 0) {
        if (ks.l == 0) continue;
        if (ks.s[0] == '#') {
            if (ks.l > 1 && ks.s[1] == '#') continue;
            header = split_ws(std::string(ks.s, ks.l));
            if (!header.empty() && header[0] == "#CHROM") {
                found = true;
                break;
            }
        } else {
            std::cerr << "wrong Line : " << std::string(ks.s, ks.l) << std::endl;
            std::cerr << "VCF Header something wrong, can't find sample info before site info" << std::endl;
            std::cerr << "VCF Header sample info Flag : [  #CHROM  ] " << std::endl;
            hts_close(fp);
            return false;
        }
    }
    hts_close(fp);
    if (!found) {
        std::cerr << "VCF Header something wrong, can't find sample info before site info" << std::endl;
        std::cerr << "VCF Header sample info Flag : [  #CHROM  ] " << std::endl;
        return false;
    }

    return resolve_subpop_header(header, specs, "[found in VCF]", "in the VCF Header", 9, out);
}

}  // namespace pld2