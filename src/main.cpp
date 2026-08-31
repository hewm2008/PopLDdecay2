#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <htslib/bgzf.h>

#include "cli.h"
#include "ehh.h"
#include "genotype_reader.h"
#include "ldkernel.h"
#include "output.h"
#include "subpop.h"
#include "types.h"
#include "vcf_reader.h"

namespace pld2 {
namespace {

std::string group_prefix(const std::string& prefix, const std::string& label) {
    return label.empty() ? prefix : prefix + "." + label;
}

}  // namespace
}  // namespace pld2

int main(int argc, char** argv) {
    using namespace pld2;
    using Clock = std::chrono::steady_clock;
    const bool timing = getenv("PLD2_TIMING") != nullptr;
    const auto ms = [](const Clock::time_point& t) {
        return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
    };

    Options opt;
    if (!parse_cli(argc, argv, opt)) return 1;

    const int max_dist_bp = opt.max_dist_kb * 1000;
    const std::string prefix = strip_stat_prefix(opt.out_stat);
    const std::string stat_path = prefix + ".stat.gz";

    // -SubPop group lists -> specs (repeated -SubPop flags and/or a two-column
    // `sample<TAB>groupid` file; a single single-column file stays legacy).
    std::vector<SubGroupSpec> specs;
    if (!parse_subpop_lists(opt.subpop_files, specs)) return 1;
    const size_t nbufs = specs.empty() ? 1 : specs.size();
    const auto label_of = [&](size_t g) {
        return specs.empty() ? std::string() : specs[g].label;
    };

    // Verify each output is writable up-front (mirrors original main()).
    for (size_t g = 0; g < nbufs; ++g) {
        const std::string p = group_prefix(prefix, label_of(g)) + ".stat.gz";
        BGZF* test = bgzf_open(p.c_str(), "w");
        if (!test) {
            std::cerr << "open OUT File error: " << p << std::endl;
            return 1;
        }
        bgzf_close(test);
    }

    ReadStats stats;
    std::vector<AccumSpec> specs_out(nbufs);
    for (AccumSpec& s : specs_out) s.out_type = opt.out_type;
    bool used_fast = false;

    // -SubPop: each group is resolved against the VCF header once (shared by
    // both paths). For a genotype-only run there is no VCF: the genotype
    // reader resolves the specs against its own #CHROM header instead, exactly
    // like the original Read_SubPopGenotype_IN.
    std::vector<SubPop> subs;
    if (!opt.in_vcf.empty() && !resolve_subpops(opt.in_vcf, specs, subs)) {
        return 1;
    }

    // Reads the -InVCF (if any) and -InGenotype (if any) inputs into packed
    // per-contig arrays (one set per group), VCF first then genotype merged
    // after, replicating the original single SNPList map per group. Returns
    // false on error.
    auto load_inputs = [&](std::vector<std::vector<ContigData>>& contigs_all,
                           ReadStats& st) -> bool {
        contigs_all.clear();
        if (!opt.in_vcf.empty()) {
            if (!read_vcf(opt.in_vcf, opt, contigs_all, st, subs)) return false;
        } else {
            contigs_all.resize(nbufs);
        }
        if (!opt.in_genotype.empty()) {
            std::vector<std::vector<GenoBuf>> gbufs_all;
            if (!read_genotype(opt.in_genotype, opt, prefix, specs, gbufs_all))
                return false;
            for (size_t g = 0; g < nbufs; ++g) {
                merge_geno_buffers(gbufs_all[g], contigs_all[g]);
            }
        }
        return true;
    };

    const auto t_run = Clock::now();

    // -EHH mode: replicate EHH_Region_LDDecay. The MAF floor of 0.05 is applied
    // before reading (mirrors the original help01), the single-pass fallback
    // reader is used (we need the per-site haplotypes of the whole file), the
    // LD/stat computation is skipped, and main returns 0 regardless of the EHH
    // run's success -- exactly like the original dispatch.
    if (!opt.ehh.empty()) {
        if (nbufs > 1) {
            std::cerr << "-EHH supports a single SubPop group only" << std::endl;
            return 1;
        }
        if (opt.maf < 0.05) {
            opt.maf = 0.05;
            std::cerr << "\t\t -MAF for -EHH should >=0.05 ; so -MAF  default 0.05 Now"
                      << std::endl;
        }
        EhhParams ep;
        if (!parse_ehh(opt.ehh, ep)) {
            std::cerr << "\tPara [-EHH] should be [chr:Site],such [chr1:5000]"
                      << std::endl;
            return 1;
        }
        std::vector<std::vector<ContigData>> contigs_all;
        stats = ReadStats{};
        if (!load_inputs(contigs_all, stats)) return 1;
        if (stats.skip_indel != 0) {
            std::cout << "warning skip Indel site, there are total skip Indel sites number is : "
                      << stats.skip_indel << std::endl;
        }
        if (stats.skip_non_biallelic != 0) {
            std::cout << "Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"
                      << stats.skip_non_biallelic << std::endl;
        }
        run_ehh(contigs_all[0], ep, max_dist_bp, opt.out_type, prefix);
        if (timing) std::cerr << "[timing] total=" << ms(t_run) << "ms\n";
        return 0;
    }

    // -OutFilterSNP forces the single-pass fallback: the filter file must list
    // passing sites in original file order (the original tool writes it during
    // its one sequential read; SubPop+OutFilterSNP writes nothing there too).
    // Genotype input also forces the fallback (no index / no fast path).
    if (opt.threads > 1 && !opt.out_filter_snp && !opt.in_vcf.empty() &&
        opt.in_genotype.empty()) {
        used_fast = compute_bins_indexed(opt.in_vcf, opt, max_dist_bp, specs_out,
                                         stats, subs);
    }
    if (!used_fast) {
        std::vector<std::vector<ContigData>> contigs_all;
        stats = ReadStats{};
        specs_out.clear();
        specs_out.resize(nbufs);
        for (AccumSpec& s : specs_out) s.out_type = opt.out_type;
        const auto t_read = Clock::now();
        if (!load_inputs(contigs_all, stats)) return 1;
        for (size_t g = 0; g < nbufs; ++g) {
            if (!contigs_all[g].empty()) {
                specs_out[g].nhaplo = contigs_all[g][0].n_haplotypes;
            }
            compute_bins(contigs_all[g], max_dist_bp, opt.threads, specs_out[g]);
        }
        const double read_ms = ms(t_read);
        const auto t_kernel = Clock::now();
        if (timing) {
            std::cerr << "[timing] fallback read+parse=" << read_ms
                      << "ms kernel=" << ms(t_kernel) << "ms\n";
        }
    }
    if (timing) std::cerr << "[timing] total=" << ms(t_run) << "ms\n";

    if (stats.skip_indel != 0) {
        std::cout << "warning skip Indel site, there are total skip Indel sites number is : "
                  << stats.skip_indel << std::endl;
    }
    if (stats.skip_non_biallelic != 0) {
        std::cout << "Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"
                  << stats.skip_non_biallelic << std::endl;
    }

    for (size_t g = 0; g < nbufs; ++g) {
        const std::string gp = group_prefix(prefix, label_of(g));
        const std::string sp = gp + ".stat.gz";
        if (opt.out_type == 4 || opt.out_type == 5) {
            if (!write_hist_stat(sp, max_dist_bp, specs_out[g].hist_rr,
                                 specs_out[g].hist_d, specs_out[g].nhaplo,
                                 opt.out_type))
                return 1;
        } else {
            if (!write_stat(sp, specs_out[g].bins, max_dist_bp, opt.out_type))
                return 1;
        }
        if (specs_out[g].ld) {
            if (!specs_out[g].ld->finish(gp)) return 1;
        }
    }

    std::cout << "Used [perl  ../bin/Plot_XX.pl ] to Plot the LDdecay" << std::endl;
    return 0;
}