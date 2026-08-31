# FORMAT — Input & Output Formats

## Input
- **VCF 4.x**: `-InVCF`, supports gzip / bgzip / plain text; BCF is also read transparently through htslib.
- Requires the GT field (`bcf_get_genotypes`); missing GT or phasing mixes follow legacy semantics.
- Files written after site filtering (`.vcf.filter.gz` from `-OutFilterSNP`) can be fed back directly into `-InVCF`.
- **InGenotype (M6)**: `-InGenotype`. **Native format** (one IUPAC token per sample, in the same style as VCF2Genotype `-WithHeader -NoRef` output):
  ```
  chr1	2310	M A A - A A A A ...
  chr1	2375	R R G R R G - R R ...
  ```
  - After `chr<TAB>site<TAB>`, one token per sample (space separated): `A/C/G/T` homozygous, `M/K/Y/R/W/S` heterozygous, `-`/`N` missing.
  - `bin/mis/plink2genotype.pl` outputs **2** allele tokens per sample and is **incompatible** with this format (do not use it for `-InGenotype`).
  - Supports gzip / plain text input.
  - For `-SubPop`, a header `#CHROM<TAB>POS<TAB>S1 S2 ...` is required (sample columns space separated, any order, matched by name).
  - Full-sample path: after FilterGeno filtering (miss/het/biallelic/MAF, Cut3base always on) **missing → heterozygous pseudo-haplotype**; SubPop path: miss/het/always-biallelic/`SeD/(SeD+Max)` MAF, **missing → (2,2)**. The two paths differ in filtering and encoding semantics (consistent with the legacy version).
- **`-InVCF` and `-InGenotype` can be given together**: sites are merged into the same per-chr table (stably sorted within the same chr, deduplicated keeping the first).

## `-L` region (M3)
- BED: `chr<TAB>start0<TAB>end1` (0-based, half-open `[start,end)`), may be `.gz`.
- Single entry: `chr` or `chr:start-end` (1-based closed interval, e.g. `chr1:1000-20000`).
- Semantics: only SNPs falling inside a region are retained; **LD is computed only between SNP pairs belonging to the same region**; multiple regions are independent, and statistics are aggregated into the same output.

## Output: `.stat.gz` (M1/M4, byte-compatible with the legacy version)
```
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
1	0.0020	NA	0.0020	NA	4
2	0.0030	NA	0.0090	NA	3
...
```
- `Dist` = physical distance (bp), 1..MaxDist(kb)*1000.
- Only bins with `count>0` are output. `Mean=Sum/count`, fixed 4 decimal places (`%.4f`).
- If a row's NumberPairs is 0, the whole row is omitted.
- With `-OutType 2/3/6/7/8` the D' column holds real values (TF=3 is always `0.0000`, since sumD is not accumulated).
- Legacy `Plot_OnePop.pl / Plot_MutiPop.pl` read it without modification.

## Output: `.stat.gz` (M4, OutType 4/5 histogram)
```
#Dist	R^2	R^2_count
#2SampleSize	120
1	0.00	3
...
```
- OT4 adds two columns `D'\tD_count` (header `#Dist\tR^2\tR^2_count\tD'\tD_count`).
- `R^2`/`D'` columns are bucket representative values `j/100` (`setprecision(2)`); `count` is the number of hits with `int(val*100)`; Dist starts at 1 and only non-zero buckets are emitted.
- `#2SampleSize` = 2× the effective sample number (per subset for SubPop).

## Output: `.vcf.filter.gz` (M3, `-OutFilterSNP`)
- Subset of VCF for sites passing the filters (Miss/Het/biallelic/MAF), format identical to the legacy version.

## Output: `.genotype.filter.gz` (M6, `-InGenotype -OutFilterSNP`)
- Subset of original genotype lines passing FilterGeno filtering (Miss/Het/biallelic/MAF + Cut3base); **plain gzip** (raw bytes == legacy ogzstream output, comparable with `cmp`); filename = `strip_stat_prefix(out)+".genotype.filter.gz"`.
- Only produced on the full-sample path (not written on the SubPop path, same as legacy); no header line.

## Output: `.LD.gz` (M4, OutType 3/6/7/8)
- Per-chr SNP×SNP pairwise LD table (plain gzip; filename = `strip_stat_prefix(out)+".LD.gz"`; row order chr lexicographic → Site1 ascending → Site2 ascending), columns:
  - OT3: `#chr	Site1	Site2	r^2	Dist`
  - OT6: `#chr	Site1	Site2	D'	r^2	Dist`
  - OT7: `#chr	Site1	Site2	D'	LOD	r^2	Dist`
  - OT8: `#chr	Site1	Site2	D'	LOD	r^2	lowCI	highCI	Dist`
- `lowCI/highCI` = 95% likelihood interval (`low_i/high_i`×1%, 101-point profile).

## Output: `.ehh.gz` (M5, `-EHH <chr:Site>`)
- EHH decay over a region centered on the core site `Site`, radius = `-MaxDist` (kb, `{prefix}` = `strip_stat_prefix(out)`); **plain gzip** (raw bytes == legacy ogzstream output, comparable with `cmp`). Triggering `-EHH` skips writing `.stat.gz`, and the success path also deletes `{prefix}.stat.gz`.
- First two header lines:
  ```
  #Chr	Site	Dist	EHH_all	EHH_0	EHH_1
  #SNP_Number	{count}
  ```
  `count` = number of non-missing sites in the region (any site with a missing sample is removed entirely).
- Data rows sorted by pos ascending over non-missing sites in the region: `chr\tSite\tDis\t{D:.4f}\t{RR:.4f}\t{(D-RR):.4f}`; sites not reaching the `D<0.088` termination threshold are zero-valued rows `…\t0.0000\t0.0000\t0.0000`.
  - `Dis = Site - coreSite` (can be negative).
  - `EHH_all = D`, `EHH_0 = RR` (contribution of the group whose first char is `'0'`), `EHH_1 = D-RR`.
- By-products: `{prefix}.ehh.pdf/.ehh.png` (plotted by Rscript, probed via `which Rscript`, only warns if missing) + `{prefix}.tmp.r` (kept at `-OutType 0`, otherwise deleted).