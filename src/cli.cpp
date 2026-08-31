#include "cli.h"

#include <cstdlib>
#include <iostream>

namespace pld2 {

namespace {

void print_short_help() {
    std::cout
        << "\n"
        << "\tUsage: PopLDdecay2 -InVCF  <in.vcf.gz>  -OutStat <out.stat>\n"
        << "\n"
        << "\t\t-InVCF         <str>     Input SNP VCF Format\n"
        << "\t\t-InGenotype    <str>     Input SNP Genotype Format\n"
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
        //<< "\t\t-Method        <int>     Method 1 (only 1 supported) [1]\n"
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
            return true;
        };

        if (flag == "InVCF" || flag == "i") {
            if (!need_arg(flag)) return false;
            opt.in_vcf = argv[++i];
        } else if (flag == "InGenotype") {
            if (!need_arg(flag)) return false;
            opt.in_genotype = argv[++i];
        } else if (flag == "SubPop" || flag == "s") {
            if (!need_arg(flag)) return false;
            opt.subpop_files.push_back(argv[++i]);
        } else if (flag == "OutStat" || flag == "o") {
            if (!need_arg(flag)) return false;
            opt.out_stat = argv[++i];
        } else if (flag == "Het") {
            if (!need_arg(flag)) return false;
            opt.het = atof(argv[++i]);
        } else if (flag == "MAF") {
            if (!need_arg(flag)) return false;
            opt.maf = atof(argv[++i]);
        } else if (flag == "Miss") {
            if (!need_arg(flag)) return false;
            opt.miss = atof(argv[++i]);
        } else if (flag == "MaxDist") {
            if (!need_arg(flag)) return false;
            opt.max_dist_kb = atoi(argv[++i]);
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
            opt.out_type = atoi(argv[++i]);
        } else if (flag == "Method") {
            if (!need_arg(flag)) return false;
            opt.method = atoi(argv[++i]);
        } else if (flag == "OutFilterSNP") {
            opt.out_filter_snp = true;
        } else if (flag == "L" || flag == "Region") {
            if (!need_arg(flag)) return false;
            opt.region = argv[++i];
            std::string err;
            if (!load_regions(opt.region, opt.regions, err)) {
                std::cerr << err << std::endl;
                return false;
            }
        } else if (flag == "T" || flag == "Threads") {
            if (!need_arg(flag)) return false;
            opt.threads = atoi(argv[++i]);
            if (opt.threads < 1) {
                std::cerr << "\t\t-T should be >= 1, reset to 1" << std::endl;
                opt.threads = 1;
            } else if (opt.threads > 256) {
                // Sanity cap: worker memory (fast-path ChrBufs, OT4/5
                // histograms) scales with thread count; clamp instead of
                // letting an accidental huge -T exhaust RAM.
                std::cerr << "\t\t-T " << opt.threads << " too large, capped to 256"
                          << std::endl;
                opt.threads = 256;
            }
        } else if (flag == "help" || flag == "h") {
            print_more_help();
            return false;
        } else {
            std::cerr << "UnKnow argument -" << flag << std::endl;
            return false;
        }
    }

    if (opt.out_stat.empty() ||
        (opt.in_vcf.empty() && opt.in_genotype.empty())) {
        std::cerr << "lack argument for the must" << std::endl;
        return false;
    }

    if (opt.out_type > 8 || opt.out_type < 0) {
        std::cerr << "\t\t-OutType should be [0-8]" << std::endl;
        return false;
    }

    if (opt.method != 1) {
        std::cerr << "\t\t-Method should be [1]; Method 2 not implemented yet" << std::endl;
        return false;
    }

    return true;
}

}  // namespace pld2