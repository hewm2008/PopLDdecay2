#include "cli.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace pld2 {

namespace {

// Only this file needs it (the integer twin is shared with ehh.cpp).
bool parse_double_strict(const char* s, double& out) {
    if (s == nullptr || *s == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (errno == ERANGE || end == s || *end != '\0') return false;
    if (!std::isfinite(v)) return false;  // "nan" / "inf" are not filter values
    out = v;
    return true;
}

// The option names this tool understands. A value that is one of them is a
// missing argument ("-OutStat -MaxDist 300" must not take "-MaxDist" as the
// output prefix); anything else, including a path that starts with '-', is
// taken as a value.
bool is_known_flag(const std::string& flag) {
    static const char* kFlags[] = {
        "InVCF",   "i",          "InGenotype", "InPlink",   "SubPop",
        "s",       "OutStat",    "o",          "Het",       "MAF",
        "Miss",    "MaxDist",    "EHH",        "OutPairLD", "OutType",
        "Method",  "OutFilterSNP", "L",        "Region",    "T",
        "Threads", "help",       "h"};
    for (const char* f : kFlags) {
        if (flag == f) return true;
    }
    return false;
}

void print_short_help() {
    std::cout
        << "\n"
        << "\tUsage: PopLDdecay2 -InVCF  <in.vcf.gz>  -OutStat <out.stat>\n"
        << "\n"
        << "\t\t-InVCF         <str>     Input SNP VCF Format\n"
        << "\t\t-InGenotype    <str>     Input SNP Genotype Format\n"
        << "\t\t-InPlink       <str>     Input PLINK [bed+bim+fam] file prefix\n"
	<< "\t\t-OutStat       <str>     OutPut Stat Dist ~ r^2/D' File\n"
        << "\n"
        << "\t\t-SubPop        <str>     SubGroup Sample File List (repeatable)\n"
        << "\t\t-MaxDist       <int>     Max Distance (kb) between two SNP [300]\n"
        << "\t\t-MAF           <float>   Min minor allele frequency filter [0.005]\n"
        << "\t\t-Het           <float>   Max ratio of het allele filter [0.88]\n"
        << "\t\t-Miss          <float>   Max ratio of miss allele filter [0.25]\n"
        << "\t\t-EHH           <str>     To Run EHH Region decay set StartSite [NA]\n"
        << "\t\t-OutFilterSNP            OutPut the final SNP to calculate\n"
        << "\t\t-L             <str>     Region restrict: BED/chr:start-end/chr\n"
        << "\t\t-T             <int>     Threads for parallel LD calculation [4]\n"
        << "\t\t-OutType       <int>     OutType [0-8], see -help\n"
        << "\t\t-help                    Show more help [hewm2008]\n"
        << "\n";
}

void print_more_help() {
    print_short_help();
    std::cout
        << "\n"
        << " Para [-OutType] can be [0-8]\n"
        << "        [-OutType 1] is the fastest for only cal (Dist ~ R^2) for MeanBin method plot\n"
        << "        [-OutType 2] will OutPut the Stat (Dist ~ r^2 & D') result for R^2 & D' MeanBin method plot\n"
        << "        [-OutType 3] will OutPut one more result of PairWise LD compaire result(with Dist~r^2)\n"
        << "        [-OutType 4] will OutPut the Stat (Dist ~ r^2 & D' ~ Number) result for R^2 & D' MeanBin/HW/MedianBin/PercentileBin plot\n"
        << "        [-OutType 5] will OutPut the Stat (Dist ~ r^2 ~ Number) result for R^2 MeanBin/HW/MedianBin/PercentileBin plot\n"
        << "        [-OutType 6] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D')\n"
        << "        [-OutType 7] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D'/LOD)\n"
        << "        [-OutType 8] will OutPut one more result of PairWise LD compaire result(with Dist~r^2/D'/LOD/CIlow/CIhi)\n"
        << "\n"
        << " Para [-EHH] format should be chr:site, such like -EHH chr1:5000000 will give out the EHH Decay of this site nearby distance\n"
        << "\n"
        << " PopLDdecay2 is a byte-compatible rewrite of PopLDdecay (hewm2008)\n"
        << " Contact:hewm2008@gmail.com / hewm2008@qq.com     join the QQ Group : 125293663\n"
        << "\n";
}

// Removes every '-' from the flag token (replicates original replace_all(flag,"-","")).
std::string strip_dashes(const char* tok) {
    std::string s(tok);
    std::string out;
    out.reserve(s.size());
    for (char c : s)
        if (c != '-') out.push_back(c);
    return out;
}

}  // namespace

bool parse_int_strict(const char* s, long long& out) {
    if (s == nullptr || *s == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const long long v = std::strtoll(s, &end, 10);
    if (errno == ERANGE || end == s || *end != '\0') return false;
    out = v;
    return true;
}

std::string strip_stat_prefix(const std::string& out_stat) {
    std::string s = out_stat;
    size_t dot = s.rfind('.');
    std::string ext = dot == std::string::npos ? s : s.substr(dot + 1);
    if (ext == "gz") {
        s = s.substr(0, s.length() - 3);
    }
    std::string base = s;
    size_t slash = s.rfind('/');
    if (slash != std::string::npos) base = s.substr(slash + 1);
    if (base != "stat") {
        dot = s.rfind('.');
        ext = dot == std::string::npos ? s : s.substr(dot + 1);
        if (ext == "stat") {
            s = s.substr(0, s.length() - 5);
        }
    }
    return s;
}

bool parse_cli(int argc, char** argv, Options& opt) {
    if (argc < 2) {
        print_short_help();
        return false;
    }

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            std::cerr << "command option error! please check." << std::endl;
            return false;
        }
        std::string flag = strip_dashes(argv[i]);

        auto need_arg = [&](const std::string& f) {
            if (i + 1 == argc) {
                std::cerr << "\t\tLack Argument for [ -" << f << " ]" << std::endl;
                return false;
            }
            // The next token is another option, not a value: taking it as one
            // would silently swallow the option (e.g. "-OutStat -MaxDist 300"
            // would name the output file "-MaxDist"). A value that merely
            // starts with '-' (a path) is still accepted.
            if (is_known_flag(strip_dashes(argv[i + 1]))) {
                std::cerr << "\t\tLack Argument for [ -" << f << " ]" << std::endl;
                return false;
            }
            return true;
        };

        if (flag == "InVCF" || flag == "i") {
            if (!need_arg(flag)) return false;
            opt.in_vcf = argv[++i];
        } else if (flag == "InGenotype") {
            if (!need_arg(flag)) return false;
            opt.in_genotype = argv[++i];
        } else if (flag == "InPlink") {
            if (!need_arg(flag)) return false;
            opt.in_plink = argv[++i];
        } else if (flag == "SubPop" || flag == "s") {
            if (!need_arg(flag)) return false;
            opt.subpop_files.push_back(argv[++i]);
        } else if (flag == "OutStat" || flag == "o") {
            if (!need_arg(flag)) return false;
            opt.out_stat = argv[++i];
        } else if (flag == "Het" || flag == "MAF" || flag == "Miss") {
            if (!need_arg(flag)) return false;
            const char* raw = argv[++i];
            double v = 0.0;
            if (!parse_double_strict(raw, v) || v < 0.0 || v > 1.0) {
                std::cerr << "\t\t-" << flag << " should be a number in [0,1], got '"
                          << raw << "'" << std::endl;
                return false;
            }
            if (flag == "Het") opt.het = v;
            else if (flag == "MAF") opt.maf = v;
            else opt.miss = v;
        } else if (flag == "MaxDist") {
            if (!need_arg(flag)) return false;
            const char* raw = argv[++i];
            long long v = 0;
            if (!parse_int_strict(raw, v)) {
                std::cerr << "\t\t-MaxDist should be an integer (kb), got '" << raw
                          << "'" << std::endl;
                return false;
            }
            if (v > 0x7FFFFFFF) v = 0x7FFFFFFF;
            opt.max_dist_kb = static_cast<int>(v);
            if (opt.max_dist_kb < 1) {
                std::cerr << "\t\t-MaxDist should be >= 1 (kb), reset to 1" << std::endl;
                opt.max_dist_kb = 1;
            } else {
                // Cap at 1 Gb minus 1 kb (1073740800 bp) in kb units. This keeps
                // max_dist_kb*1000 well below INT32_MAX (distance bins are
                // int64-indexed, but max_dist_bp and the OT4/5 histogram slot
                // math are int-sized), mirroring the -T clamp: an accidental
                // huge -MaxDist must not overflow / exhaust RAM.
                const int kMaxDistKbCap = ((1 << 30) - 1024) / 1000;
                if (opt.max_dist_kb > kMaxDistKbCap) {
                    std::cerr << "\t\t-MaxDist " << opt.max_dist_kb
                              << " too large, capped to " << kMaxDistKbCap
                              << " kb (1 Gb - 1 kb)" << std::endl;
                    opt.max_dist_kb = kMaxDistKbCap;
                }
            }
        } else if (flag == "EHH") {
            if (!need_arg(flag)) return false;
            opt.ehh = argv[++i];
            if (opt.ehh.find(':') == std::string::npos) {
                std::cerr << "\tPara [-EHH] should be [chr:Site],such [chr1:5000]"
                          << std::endl;
                return false;
            }
        } else if (flag == "OutPairLD" || flag == "OutType") {
            if (!need_arg(flag)) return false;
            const char* raw = argv[++i];
            long long v = 0;
            if (!parse_int_strict(raw, v) || v < 0 || v > 8) {
                std::cerr << "\t\t-OutType should be an integer in [0-8], got '"
                          << raw << "'" << std::endl;
                return false;
            }
            opt.out_type = static_cast<int>(v);
        } else if (flag == "Method") {
            if (!need_arg(flag)) return false;
            const char* raw = argv[++i];
            long long v = 0;
            if (!parse_int_strict(raw, v)) {
                std::cerr << "\t\t-Method should be an integer, got '" << raw << "'"
                          << std::endl;
                return false;
            }
            opt.method = static_cast<int>(v);
        } else if (flag == "OutFilterSNP") {
            opt.out_filter_snp = true;
        } else if (flag == "L" || flag == "Region") {
            if (!need_arg(flag)) return false;
            const std::string region = argv[++i];
            std::string err;
            if (!load_regions(region, opt.regions, err)) {
                std::cerr << err << std::endl;
                return false;
            }
        } else if (flag == "T" || flag == "Threads") {
            if (!need_arg(flag)) return false;
            const char* raw = argv[++i];
            long long v = 0;
            if (!parse_int_strict(raw, v)) {
                std::cerr << "\t\t-T should be an integer, got '" << raw << "'"
                          << std::endl;
                return false;
            }
            if (v < 1) {
                // Sanity: worker memory (fast-path ChrBufs, OT4/5 histograms)
                // scales with the thread count.
                std::cerr << "\t\t-T should be >= 1, reset to 1" << std::endl;
                v = 1;
            } else if (v > 256) {
                std::cerr << "\t\t-T " << v << " too large, capped to 256"
                          << std::endl;
                v = 256;
            }
            opt.threads = static_cast<int>(v);
        } else if (flag == "help" || flag == "h") {
            print_more_help();
            return false;
        } else {
            std::cerr << "UnKnow argument -" << flag << std::endl;
            return false;
        }
    }

    if (opt.out_stat.empty() ||
        (opt.in_vcf.empty() && opt.in_genotype.empty() && opt.in_plink.empty())) {
        std::cerr << "lack argument for the must" << std::endl;
        return false;
    }

    // PLINK carries its own sample list and genotype coding, so it cannot be
    // merged with the other two inputs (their sample order and set cannot be
    // cross-checked). -InVCF + -InGenotype remain combinable, as before.
    if (!opt.in_plink.empty() && (!opt.in_vcf.empty() || !opt.in_genotype.empty())) {
        std::cerr << "\t\t[-InPlink] can't be used with [-InVCF] or [-InGenotype]"
                  << std::endl;
        return false;
    }

    if (opt.method != 1) {
        std::cerr << "\t\t-Method should be [1]; Method 2 not implemented yet" << std::endl;
        return false;
    }

    return true;
}

}  // namespace pld2
