#include "ldkernel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>

namespace pld2 {

namespace {

// LN10 = log(10.0), computed once exactly like the original statementVar ctor.
const double kLN10 = std::log(10.0);

// ---------------------------------------------------------------------------
// Per-pair formula replicas. Each reproduces the original cal_RR_* op-for-op
// (Calculate.h) so double rounding is identical -> byte-identical output.
// The 2x2 haplotype table is (n00=dde[0][0], n01=dde[0][1], n10=dde[1][0],
// n11=dde[1][1]); ALL_count in the original is a double.
// ---------------------------------------------------------------------------

// cal_RR_MA (used by the original for OutType 3 and 5): parenthesized
// Cal_A and r^2 = (D_A*D_A)/(Cal_A*Cal_B).
inline bool cal_rr_ma(int n11, int n10, int n01, int n00, double& rr) {
    const int tmpAA = n11 + n10;
    if (tmpAA == 0) return false;
    if ((n11 + n01) == 0) return false;
    const double all_count = static_cast<double>(n00 + n01 + n10 + n11);
    double prob0 = static_cast<double>(n00) / all_count;
    double prob1 = static_cast<double>(n01) / all_count;
    double prob2 = static_cast<double>(n10) / all_count;
    double pA1 = prob0 + prob1;
    double pA2 = prob0 + prob2;
    double cal_b = pA1 * pA2;
    double cal_a = 1.0 - (pA1 + pA2) + cal_b;
    if (cal_a == 0.0 || cal_b == 0.0) {
        if (prob0 < 1e-10) prob0 = 1e-10;
        if (prob1 < 1e-10) prob1 = 1e-10;
        if (prob2 < 1e-10) prob2 = 1e-10;
        pA1 = prob0 + prob1;
        pA2 = prob0 + prob2;
        cal_b = pA1 * pA2;
        cal_a = 1.0 - (pA1 + pA2) + cal_b;
    }
    const double d_a = prob0 - cal_b;
    rr = (d_a * d_a) / (cal_a * cal_b);
    return true;
}

// Inline OutType-1 r^2 (the original ProMethod1.h TF==1 block). Differs from
// cal_RR_MA in exactly one spot: Cal_A uses the un-parenthesized
// 1.0-pA1-pA2+Cal_B (r^2 itself is the same (D_A*D_A)/(Cal_A*Cal_B)).
inline bool cal_rr_ot1(int n11, int n10, int n01, int n00, double& rr) {
    const int tmpAA = n11 + n10;
    if (tmpAA == 0) return false;
    if ((n11 + n01) == 0) return false;
    const double all_count = static_cast<double>(n00 + n01 + n10 + n11);
    double prob0 = static_cast<double>(n00) / all_count;
    double prob1 = static_cast<double>(n01) / all_count;
    double prob2 = static_cast<double>(n10) / all_count;
    double pA1 = prob0 + prob1;
    double pA2 = prob0 + prob2;
    double cal_b = pA1 * pA2;
    double cal_a = 1.0 - pA1 - pA2 + cal_b;
    if (cal_a == 0.0 || cal_b == 0.0) {
        if (prob0 < 1e-10) prob0 = 1e-10;
        if (prob1 < 1e-10) prob1 = 1e-10;
        if (prob2 < 1e-10) prob2 = 1e-10;
        pA1 = prob0 + prob1;
        pA2 = prob0 + prob2;
        cal_b = pA1 * pA2;
        cal_a = 1.0 - pA1 - pA2 + cal_b;
    }
    const double d_a = prob0 - cal_b;
    rr = (d_a * d_a) / (cal_a * cal_b);
    return true;
}

// cal_RR_D_MA.
inline bool cal_rr_d_ma(int n11, int n10, int n01, int n00, double& d,
                        double& rr) {
    const int tmpAA = n11 + n10;
    if (tmpAA == 0) return false;
    if ((n11 + n01) == 0) return false;
    const double all_count = static_cast<double>(n00 + n01 + n10 + n11);
    double prob0 = static_cast<double>(n00) / all_count;
    double prob1 = static_cast<double>(n01) / all_count;
    double prob2 = static_cast<double>(n10) / all_count;
    double pA1 = prob0 + prob1;
    double pB1 = 1.0 - pA1;
    double pA2 = prob0 + prob2;
    double pB2 = 1.0 - pA2;
    double xpA1_pA2 = pA1 * pA2;
    double d_a = prob0 - xpA1_pA2;
    double cal_a, cal_b;
    if (d_a > 0) {
        cal_a = pB1 * pA2;
        cal_b = pA1 * pB2;
    } else {
        d_a = 0.0 - d_a;
        cal_a = pB1 * pB2;
        cal_b = xpA1_pA2;
    }
    double d_max = cal_a;
    if (cal_a > cal_b) d_max = cal_b;
    if (d_max == 0) {
        if (prob0 < 1e-10) prob0 = 1e-10;
        if (prob1 < 1e-10) prob1 = 1e-10;
        if (prob2 < 1e-10) prob2 = 1e-10;
        pA1 = prob0 + prob1;
        pB1 = 1.0 - pA1;
        pA2 = prob0 + prob2;
        pB2 = 1.0 - pA2;
        xpA1_pA2 = pA1 * pA2;
        d_a = prob0 - xpA1_pA2;
        if (d_a > 0) {
            cal_a = pB1 * pA2;
            cal_b = pA1 * pB2;
        } else {
            d_a = 0.0 - d_a;
            cal_a = pB1 * pB2;
            cal_b = xpA1_pA2;
        }
        d_max = cal_a;
        if (cal_a > cal_b) d_max = cal_b;
    }
    d = d_a / d_max;
    rr = (d_a / cal_a) * (d_a / cal_b);
    return true;
}

// cal_RR_D2_MA.
inline bool cal_rr_d2_ma(int n11, int n10, int n01, int n00, double& d,
                         double& rr, double& lod) {
    double known[4];
    known[0] = static_cast<double>(n00);
    known[1] = static_cast<double>(n01);
    known[2] = static_cast<double>(n10);
    known[3] = static_cast<double>(n11);
    const double tmpAA = known[3] + known[2];
    if (tmpAA == 0) return false;
    if ((known[3] + known[1]) == 0) return false;
    const double all_count = known[0] + known[1] + known[2] + known[3];
    double prob0 = known[0] / all_count;
    double prob1 = known[1] / all_count;
    double prob2 = known[2] / all_count;
    double prob3 = 1.0 - prob0 - prob1 - prob2;
    if (prob0 < 1e-10) prob0 = 1e-10;
    if (prob1 < 1e-10) prob1 = 1e-10;
    if (prob2 < 1e-10) prob2 = 1e-10;
    if (prob3 < 1e-10) prob3 = 1e-10;
    const double pA1 = prob0 + prob1;
    const double pB1 = 1.0 - pA1;
    const double pA2 = prob0 + prob2;
    const double pB2 = 1.0 - pA2;
    const double xpA1_pA2 = pA1 * pA2;
    const double xpA1_pB2 = pA1 * pB2;
    const double xpB1_pA2 = pB1 * pA2;
    const double xpB1_pB2 = pB1 * pB2;
    const double loglike1 = (known[0] * std::log(prob0) +
                             known[1] * std::log(prob1) +
                             known[2] * std::log(prob2) +
                             known[3] * std::log(prob3)) /
                            kLN10;
    const double loglike0 = (known[0] * std::log(xpA1_pA2) +
                             known[1] * std::log(xpA1_pB2) +
                             known[2] * std::log(xpB1_pA2) +
                             known[3] * std::log(xpB1_pB2)) /
                            kLN10;
    double d_a = prob0 - xpA1_pA2;
    double cal_a, cal_b;
    if (d_a < 0) {
        d_a = 0.0 - d_a;
        cal_a = xpB1_pB2;
        cal_b = xpA1_pA2;
    } else {
        cal_a = xpB1_pA2;
        cal_b = xpA1_pB2;
    }
    double d_max = cal_a;
    if (cal_a > cal_b) d_max = cal_b;
    d = d_a / d_max;
    rr = (d_a / cal_a) * (d_a / cal_b);
    lod = loglike1 - loglike0;
    return true;
}

// cal_RR_D3_MA (includes the 101-point likelihood-surface CI).
inline bool cal_rr_d3_ma(int n11, int n10, int n01, int n00, double& d,
                         double& rr, double& lod, int& low_i, int& high_i) {
    double known[4];
    known[0] = static_cast<double>(n00);
    known[1] = static_cast<double>(n01);
    known[2] = static_cast<double>(n10);
    known[3] = static_cast<double>(n11);
    const double tmpAA = known[3] + known[2];
    if (tmpAA == 0) return false;
    if ((known[3] + known[1]) == 0) return false;
    const double all_count = known[0] + known[1] + known[2] + known[3];
    double prob[4];
    prob[0] = known[0] / all_count;
    prob[1] = known[1] / all_count;
    prob[2] = known[2] / all_count;
    prob[3] = 1.0 - prob[0] - prob[1] - prob[2];
    if (prob[0] < 1e-10) prob[0] = 1e-10;
    if (prob[1] < 1e-10) prob[1] = 1e-10;
    if (prob[2] < 1e-10) prob[2] = 1e-10;
    if (prob[3] < 1e-10) prob[3] = 1e-10;
    double pA1 = prob[0] + prob[1];
    const double pB1 = 1.0 - pA1;
    double pA2 = prob[0] + prob[2];
    double pB2 = 1.0 - pA2;
    double xpA1_pA2 = pA1 * pA2;
    const double xpA1_pB2 = pA1 * pB2;
    const double xpB1_pA2 = pB1 * pA2;
    const double xpB1_pB2 = pB1 * pB2;
    double loglike1 = (known[0] * std::log(prob[0]) +
                       known[1] * std::log(prob[1]) +
                       known[2] * std::log(prob[2]) +
                       known[3] * std::log(prob[3])) /
                      kLN10;
    const double loglike0 = (known[0] * std::log(xpA1_pA2) +
                             known[1] * std::log(xpA1_pB2) +
                             known[2] * std::log(xpB1_pA2) +
                             known[3] * std::log(xpB1_pB2)) /
                            kLN10;
    double d_a = prob[0] - xpA1_pA2;
    double cal_a, cal_b;
    if (d_a < 0) {
        // Swap prob(0<->1), prob(3<->2), pA2<->pB2 and known(0<->1),
        // known(3<->2); replicate using a temp like the original's ALL_count.
        double t = prob[0];
        prob[0] = prob[1];
        prob[1] = t;
        t = prob[3];
        prob[3] = prob[2];
        prob[2] = t;
        pA2 = pA2 + pB2;
        pB2 = pA2 - pB2;
        pA2 = pA2 - pB2;
        d_a = 0.0 - d_a;
        t = known[0];
        known[0] = known[1];
        known[1] = t;
        t = known[3];
        known[3] = known[2];
        known[2] = t;
        cal_a = pA2 * pB1;
        cal_b = pA1 * pB2;
    } else {
        cal_a = xpB1_pA2;
        cal_b = xpA1_pB2;
    }
    double d_max = cal_a;
    if (cal_a > cal_b) d_max = cal_b;
    d = d_a / d_max;
    rr = (d_a / cal_a) * (d_a / cal_b);
    lod = loglike1 - loglike0;

    // 101-point likelihood surface for the 5% CI.
    xpA1_pA2 = pA1 * pA2;
    double lsurface[101];
    double tmp_aa, tmp_ab, tmp_ba, tmp_bb, dpr;
    for (int i = 0; i < 100; i++) {
        dpr = static_cast<double>(i) * 0.01;
        tmp_aa = dpr * d_max + xpA1_pA2;
        tmp_ab = pA1 - tmp_aa;
        tmp_ba = pA2 - tmp_aa;
        tmp_bb = pB1 - tmp_ba;
        lsurface[i] = (known[0] * std::log(tmp_aa) + known[1] * std::log(tmp_ab) +
                       known[2] * std::log(tmp_ba) + known[3] * std::log(tmp_bb)) /
                      kLN10;
    }
    dpr = static_cast<double>(100) * 0.01;
    tmp_aa = dpr * d_max + xpA1_pA2;
    tmp_ab = pA1 - tmp_aa;
    tmp_ba = pA2 - tmp_aa;
    tmp_bb = pB1 - tmp_ba;
    if (tmp_aa < 1e-10) tmp_aa = 1e-10;
    if (tmp_ab < 1e-10) tmp_ab = 1e-10;
    if (tmp_ba < 1e-10) tmp_ba = 1e-10;
    if (tmp_bb < 1e-10) tmp_bb = 1e-10;
    lsurface[100] = (known[0] * std::log(tmp_aa) + known[1] * std::log(tmp_ab) +
                     known[2] * std::log(tmp_ba) + known[3] * std::log(tmp_bb)) /
                    kLN10;

    double total_prob = 0.0;
    for (int i = 0; i <= 100; i++) {
        lsurface[i] -= loglike1;
        lsurface[i] = std::pow(10.0, lsurface[i]);
        total_prob += lsurface[i];
    }
    const double cut5off = total_prob * 0.05;

    double sum_prob = 0.0;
    low_i = 0;
    for (int i = 0; i <= 100; i++) {
        sum_prob += lsurface[i];
        if (sum_prob > cut5off && (sum_prob - lsurface[i]) < cut5off) {
            low_i = i - 1;
            break;
        }
    }
    sum_prob = 0.0;
    high_i = 0;
    for (int i = 100; i >= 0; i--) {
        sum_prob += lsurface[i];
        if (sum_prob > cut5off && (sum_prob - lsurface[i]) < cut5off) {
            high_i = i + 1;
            break;
        }
    }
    if (high_i > 100) high_i = 100;
    return true;
}

// Bucket index for an OT4/5 value histogram: 101 slots per distance bin.
// size_t arithmetic so a large -MaxDist (up to the 1 Gb - 1 kb CLI cap) cannot
// overflow; the histogram vectors are indexed by size_t everywhere already.
inline size_t hist_slot(int64_t dis, double val) {
    return static_cast<size_t>(dis) * 101 + static_cast<size_t>(val * 100.0);
}

// Dispatches one validated pair to the accumulator(s) implied by out_type.
inline void accumulate_pair(const ContigData& cd, size_t i, size_t j,
                            int64_t dis, int n11, int n10, int n01, int n00,
                            AccumSpec& spec) {
    const int ot = spec.out_type;
    if (ot == 0) return;
    if (ot == 1) {
        double rr;
        if (!cal_rr_ot1(n11, n10, n01, n00, rr)) return;
        Bin& bin = spec.bins[static_cast<size_t>(dis)];
        bin.count++;
        bin.sumRR += rr;
        return;
    }
    if (ot == 3 || ot == 5) {
        double rr;
        if (!cal_rr_ma(n11, n10, n01, n00, rr)) return;
        if (ot == 5) {
            if (!spec.hist_rr.empty()) spec.hist_rr[hist_slot(dis, rr)]++;
            return;
        }
        Bin& bin = spec.bins[static_cast<size_t>(dis)];
        bin.count++;
        bin.sumRR += rr;
        if (ot == 3 && spec.ld) {
            spec.ld->emit_rr(cd.name, cd.snps[i].pos, cd.snps[j].pos, rr, dis);
        }
        return;
    }
    if (ot == 2 || ot == 4 || ot == 6) {
        double d, rr;
        if (!cal_rr_d_ma(n11, n10, n01, n00, d, rr)) return;
        if (ot == 4) {
            if (!spec.hist_rr.empty()) spec.hist_rr[hist_slot(dis, rr)]++;
            if (!spec.hist_d.empty()) spec.hist_d[hist_slot(dis, d)]++;
            return;
        }
        Bin& bin = spec.bins[static_cast<size_t>(dis)];
        bin.count++;
        bin.sumRR += rr;
        bin.sumD += d;
        if (ot == 6 && spec.ld) {
            spec.ld->emit_d_rr(cd.name, cd.snps[i].pos, cd.snps[j].pos, d, rr, dis);
        }
        return;
    }
    if (ot == 7) {
        double d, rr, lod;
        if (!cal_rr_d2_ma(n11, n10, n01, n00, d, rr, lod)) return;
        Bin& bin = spec.bins[static_cast<size_t>(dis)];
        bin.count++;
        bin.sumRR += rr;
        bin.sumD += d;
        if (spec.ld) {
            spec.ld->emit_d_lod_rr(cd.name, cd.snps[i].pos, cd.snps[j].pos, d, lod, rr, dis);
        }
        return;
    }
    // ot == 8
    double d, rr, lod;
    int low_i, high_i;
    if (!cal_rr_d3_ma(n11, n10, n01, n00, d, rr, lod, low_i, high_i)) return;
    Bin& bin = spec.bins[static_cast<size_t>(dis)];
    bin.count++;
    bin.sumRR += rr;
    bin.sumD += d;
    if (spec.ld) {
        spec.ld->emit_d_lod_rr_ci(cd.name, cd.snps[i].pos, cd.snps[j].pos, d, lod, rr,
                                  low_i, high_i, dis);
    }
}

// nwords==2 specialization: hoists the (i-site) bitmaps + hapmask out of the
// j-loop and unrolls the 2 words, cutting per-pair loop/load overhead ~1.4x.
// Arithmetic is identical to the generic loop -> byte-identical bins.
//
// OutType 1 (the default) is specialized via the template: cal_rr_ot1 + the
// bin update are inlined directly in the pair loop instead of calling
// accumulate_pair. (The M4 accumulate_pair refactor cost the 0.3.x inline
// loop ~45%: bb T1 kernel 101s -> 147s; profiling the binary shows a
// callq to accumulate_pair per pair.) Other OutTypes keep the shared
// dispatch path.
template <bool kOt1Inline>
void accumulate_contig_nwords2_t(const ContigData& cd, int max_dist_bp,
                                 AccumSpec& spec, size_t core_n) {
    const size_t n = cd.snps.size();
    const uint64_t m0 = ~uint64_t(0);
    const size_t last_bits = cd.n_haplotypes % 64;
    const uint64_t m1 =
        (last_bits == 0) ? ~uint64_t(0) : ((uint64_t(1) << last_bits) - 1);

    for (size_t i = 0; i < core_n; ++i) {
        const SnpData& a = cd.snps[i];
        const int64_t pos_i = a.pos;
        const uint64_t* amin = a.minor.data();
        const uint64_t* amis = a.miss.data();
        const uint64_t a0 = amin[0], a1 = amin[1];
        const uint64_t am0 = amis[0], am1 = amis[1];
        for (size_t j = i + 1; j < n; ++j) {
            const SnpData& b = cd.snps[j];
            const int64_t dis = b.pos - pos_i;
            if (dis > max_dist_bp) break;
            const uint64_t* bmin = b.minor.data();
            const uint64_t* bmis = b.miss.data();
            const uint64_t b0 = bmin[0], b1 = bmin[1];
            const uint64_t bm0 = bmis[0], bm1 = bmis[1];
            const uint64_t v0 = ~(am0 | bm0);
            const uint64_t v1 = ~(am1 | bm1);
            const int n11 = __builtin_popcountll(a0 & b0 & m0) +
                            __builtin_popcountll(a1 & b1 & m1);
            const int n10 = __builtin_popcountll(a0 & ~b0 & v0 & m0) +
                            __builtin_popcountll(a1 & ~b1 & v1 & m1);
            const int n01 = __builtin_popcountll(~a0 & b0 & v0 & m0) +
                            __builtin_popcountll(~a1 & b1 & v1 & m1);
            const int n00 = __builtin_popcountll(~a0 & ~b0 & v0 & m0) +
                            __builtin_popcountll(~a1 & ~b1 & v1 & m1);
            if constexpr (kOt1Inline) {
                double rr;
                if (!cal_rr_ot1(n11, n10, n01, n00, rr)) continue;
                Bin& bin = spec.bins[static_cast<size_t>(dis)];
                bin.count++;
                bin.sumRR += rr;
            } else {
                accumulate_pair(cd, i, j, dis, n11, n10, n01, n00, spec);
            }
        }
    }
}

void accumulate_contig_nwords2(const ContigData& cd, int max_dist_bp,
                               AccumSpec& spec, size_t core_n) {
    if (spec.out_type == 1) {
        accumulate_contig_nwords2_t<true>(cd, max_dist_bp, spec, core_n);
    } else {
        accumulate_contig_nwords2_t<false>(cd, max_dist_bp, spec, core_n);
    }
}

// nwords==1 specialization (N<=32 samples): the single bitmap word is hoisted
// out of the j-loop. Arithmetic is identical to the generic loop.
template <bool kOt1Inline>
void accumulate_contig_nwords1_t(const ContigData& cd, int max_dist_bp,
                                 AccumSpec& spec, size_t core_n) {
    const size_t n = cd.snps.size();
    const size_t last_bits = cd.n_haplotypes % 64;
    const uint64_t m0 =
        (last_bits == 0) ? ~uint64_t(0) : ((uint64_t(1) << last_bits) - 1);

    for (size_t i = 0; i < core_n; ++i) {
        const SnpData& a = cd.snps[i];
        const int64_t pos_i = a.pos;
        const uint64_t* amin = a.minor.data();
        const uint64_t* amis = a.miss.data();
        const uint64_t a0 = amin[0];
        const uint64_t am0 = amis[0];
        for (size_t j = i + 1; j < n; ++j) {
            const SnpData& b = cd.snps[j];
            const int64_t dis = b.pos - pos_i;
            if (dis > max_dist_bp) break;
            const uint64_t* bmin = b.minor.data();
            const uint64_t* bmis = b.miss.data();
            const uint64_t b0 = bmin[0];
            const uint64_t bm0 = bmis[0];
            const uint64_t v0 = ~(am0 | bm0);
            const int n11 = __builtin_popcountll(a0 & b0 & m0);
            const int n10 = __builtin_popcountll(a0 & ~b0 & v0 & m0);
            const int n01 = __builtin_popcountll(~a0 & b0 & v0 & m0);
            const int n00 = __builtin_popcountll(~a0 & ~b0 & v0 & m0);
            if constexpr (kOt1Inline) {
                double rr;
                if (!cal_rr_ot1(n11, n10, n01, n00, rr)) continue;
                Bin& bin = spec.bins[static_cast<size_t>(dis)];
                bin.count++;
                bin.sumRR += rr;
            } else {
                accumulate_pair(cd, i, j, dis, n11, n10, n01, n00, spec);
            }
        }
    }
}

void accumulate_contig_nwords1(const ContigData& cd, int max_dist_bp,
                               AccumSpec& spec, size_t core_n) {
    if (spec.out_type == 1) {
        accumulate_contig_nwords1_t<true>(cd, max_dist_bp, spec, core_n);
    } else {
        accumulate_contig_nwords1_t<false>(cd, max_dist_bp, spec, core_n);
    }
}

void accumulate_contig_impl(const ContigData& cd, int max_dist_bp,
                            AccumSpec& spec, size_t core_n) {
    const size_t n = cd.snps.size();
    if (n == 0) return;
    if (core_n > n) core_n = n;
    const size_t nwords = cd.nwords;
    const size_t last_bits = cd.n_haplotypes % 64;

    if (nwords == 1) {
        accumulate_contig_nwords1(cd, max_dist_bp, spec, core_n);
        return;
    }
    if (nwords == 2) {
        accumulate_contig_nwords2(cd, max_dist_bp, spec, core_n);
        return;
    }

    for (size_t i = 0; i < core_n; ++i) {
        const SnpData& a = cd.snps[i];
        const int64_t pos_i = a.pos;
        const uint64_t* amin = a.minor.data();
        const uint64_t* amis = a.miss.data();
        for (size_t j = i + 1; j < n; ++j) {
            const SnpData& b = cd.snps[j];
            const int64_t dis = b.pos - pos_i;
            if (dis > max_dist_bp) break;
            const uint64_t* bmin = b.minor.data();
            const uint64_t* bmis = b.miss.data();

            // 2x2 haplotype table via 4 popcounts per word.
            // n00=dde[0][0], n01=dde[0][1], n10=dde[1][0], n11=dde[1][1]
            int n11 = 0, n10 = 0, n01 = 0, n00 = 0;
            for (size_t w = 0; w < nwords; ++w) {
                const uint64_t ai = amin[w];
                const uint64_t bi = bmin[w];
                const uint64_t valid = ~(amis[w] | bmis[w]);
                // Mask off padding bits beyond the real 2*N haplotypes.
                const uint64_t hapmask =
                    (w + 1 < nwords)
                        ? ~uint64_t(0)
                        : (last_bits == 0 ? ~uint64_t(0)
                                          : ((uint64_t(1) << last_bits) - 1));
                n11 += __builtin_popcountll(ai & bi & hapmask);
                n10 += __builtin_popcountll(ai & ~bi & valid & hapmask);
                n01 += __builtin_popcountll(~ai & bi & valid & hapmask);
                n00 += __builtin_popcountll(~ai & ~bi & valid & hapmask);
            }

            accumulate_pair(cd, i, j, dis, n11, n10, n01, n00, spec);
        }
    }
}

}  // namespace

void accumulate_contig(const ContigData& cd, int max_dist_bp, AccumSpec& spec,
                       size_t core_n, const std::string& part_label) {
    if (spec.ld) {
        spec.ld->begin_contig(part_label.empty() ? cd.name : part_label);
    }
    accumulate_contig_impl(cd, max_dist_bp, spec, core_n);
}

void compute_bins(const std::vector<ContigData>& contigs, int max_dist_bp,
                  int nthreads, AccumSpec& out) {
    const size_t nbins = static_cast<size_t>(max_dist_bp) + 1;
    const bool hist_mode = (out.out_type == 4 || out.out_type == 5);
    out.bins.assign(nbins, Bin{});
    if (hist_mode) {
        out.hist_rr.assign(nbins * 101, 0);
        out.hist_d.assign((out.out_type == 4) ? nbins * 101 : 0, 0);
    }
    if (out.out_type == 3 || out.out_type == 6 || out.out_type == 7 ||
        out.out_type == 8) {
        out.ld = std::make_unique<LdWriter>(out.out_type);
    }

    std::vector<size_t> work;
    for (size_t i = 0; i < contigs.size(); ++i) {
        if (!contigs[i].snps.empty()) work.push_back(i);
    }
    if (work.empty()) return;

    // Sequential path: byte-identical accumulation order to the original.
    if (nthreads <= 1 || work.size() == 1) {
        for (size_t idx : work) accumulate_contig(contigs[idx], max_dist_bp, out);
        return;
    }

    // Greedy largest-first assignment to the least-loaded worker.
    size_t nworkers = std::min<size_t>(static_cast<size_t>(nthreads), work.size());
    if (hist_mode && nworkers > 4) {
        // OutType 4/5 keeps a (MaxDist+1)*101 int histogram per worker
        // (~121 MB/worker at MaxDist=300 kb). Cap the worker count to bound
        // memory; integer sums merge in any order, so output is byte-identical.
        std::cerr << "OutType 4/5 histogram mode: capping accumulation workers "
                     "to 4 (histogram is ~121 MB/worker at MaxDist=300 kb)\n";
        nworkers = 4;
    }
    std::vector<size_t> order = work;
    std::sort(order.begin(), order.end(), [&](size_t x, size_t y) {
        return contigs[x].snps.size() > contigs[y].snps.size();
    });
    std::vector<std::vector<size_t>> assign(nworkers);
    std::vector<uint64_t> load(nworkers, 0);
    for (size_t idx : order) {
        size_t w = static_cast<size_t>(std::min_element(load.begin(), load.end()) - load.begin());
        assign[w].push_back(idx);
        load[w] += static_cast<uint64_t>(contigs[idx].snps.size());
    }

    std::vector<AccumSpec> partial(nworkers);
    for (auto& p : partial) {
        p.out_type = out.out_type;
        p.nhaplo = out.nhaplo;
        p.bins.assign(nbins, Bin{});
        if (hist_mode) {
            p.hist_rr.assign(nbins * 101, 0);
            if (out.out_type == 4) p.hist_d.assign(nbins * 101, 0);
        }
        if (out.ld) p.ld = std::make_unique<LdWriter>(out.out_type);
    }

    std::vector<std::thread> pool;
    pool.reserve(nworkers - 1);
    for (size_t w = 1; w < nworkers; ++w) {
        pool.emplace_back([&, w] {
            for (size_t idx : assign[w]) {
                accumulate_contig(contigs[idx], max_dist_bp, partial[w]);
            }
        });
    }
    for (size_t idx : assign[0]) accumulate_contig(contigs[idx], max_dist_bp, partial[0]);
    for (auto& t : pool) t.join();

    // Merge partials in worker order (deterministic). Integer counts exact;
    // double sums ordered per worker (see ldkernel.h note).
    for (size_t w = 0; w < nworkers; ++w) {
        for (size_t d = 0; d < nbins; ++d) {
            out.bins[d].count += partial[w].bins[d].count;
            out.bins[d].sumRR += partial[w].bins[d].sumRR;
            out.bins[d].sumD += partial[w].bins[d].sumD;
        }
        if (hist_mode) {
            for (size_t k = 0; k < nbins * 101; ++k) {
                out.hist_rr[k] += partial[w].hist_rr[k];
                if (out.out_type == 4) out.hist_d[k] += partial[w].hist_d[k];
            }
        }
        if (out.ld && partial[w].ld) out.ld->absorb(*partial[w].ld);
    }
}

}  // namespace pld2
