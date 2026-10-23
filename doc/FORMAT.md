# FORMAT — input and output formats

## Command-line parameter ranges
| Parameter | Accepted value | On an out-of-range or invalid value |
|---|---|---|
| `-MAF` / `-Het` / `-Miss` | A finite decimal in `[0, 1]` (both ends included) | Error, exit code 1; the message names the option and the offending value |
| `-MaxDist` | An integer (kb) | A non-integer is an error; `< 1` is still reset to 1 with a message (as before); a huge value (> 1 Gb − 1 kb) is capped, with a message |
| `-OutType` | An integer, `0`–`8` | Out of range or not an integer: error |
| `-Method` | An integer; only `1` is implemented | Any other value: error (same wording as the original: `Method 2 not implemented yet`) |
| `-T` / `-Threads` | An integer, `1`–`256` | A non-integer is an error; `< 1` is reset to 1 and `> 256` is capped, both with a message |
| `-EHH` | `chr:integer` | A missing `:` or a non-integer site: error |

- **The whole token has to be a number.** `atoi()`/`atof()` cannot tell `"abc"`
  from `0`, so the original silently turned `-MAF abc` into 0 (MAF filter off),
  `-OutType abc` into 0 (a header-only table), `-MaxDist 1e9` into 1 kb and
  `-EHH chr1:abc` into a core site at 0. All of these are now reported.
- The boundary values stay valid: `-MAF 0 -Miss 1 -Het 1` (every filter wide
  open) and `-MaxDist 0` (reset to 1).
- **Missing argument**: if the next token is another known option name, the
  argument is missing (`-OutStat -MaxDist 300` no longer names the output file
  `-MaxDist`). A value that merely starts with `-` is still accepted, so a path
  may start with a dash.
- Relation to the original: the original uses `atoi`/`atof` too and therefore
  accepts these as well. The difference is confined to the error path (the
  message and the exit code); a successful run writes the same bytes. A failed
  `-EHH` run still exits 0, exactly like 3.45 -- a constraint that a CI case
  now pins.

## Output file overwrite semantics
- At startup the target `{prefix}.stat.gz` is opened for writing to check that
  it is writable -- **which truncates a previous result with the same name**,
  even if the run then fails on the input. PopLDdecay 3.45 does exactly the
  same (it opens the same file with an `ogzstream` for the same check), so this
  is a deliberate replication, not an oversight.
- Likewise, with populations (`-SubPop`) every group's
  `{prefix}.{label}.stat.gz` is truncated.
- `.LD.gz` / `.ehh.gz` / `.vcf.filter.gz` / `.genotype.filter.gz` /
  `.plink.filter.bim` are written only when the corresponding feature is
  enabled; a failed write leaves no partial file behind (it is removed and the
  exit code is non-zero).

## Input
- **VCF**: `-InVCF`, plain / gzip / bgzip. BCF is **not** supported: the reader
  reads the file as text lines.
- A GT field is required; a file without GT, or with mixed phasing, is handled
  the way the original handles it.
- The filtered file written by `-OutFilterSNP` (`.vcf.filter.gz`) can be fed
  straight back into `-InVCF`.
- **InGenotype (M6)**: `-InGenotype`. The **native format** (one IUPAC code
  per sample, in the style of VCF2Genotype `-WithHeader -NoRef`):
  ```
  chr1	2310	M A A - A A A A ...
  chr1	2375	R R G R R G - R R ...
  ```
  - `chr<TAB>site<TAB>` followed by one code per sample (space separated):
    `A/C/G/T` homozygous, `M/K/Y/R/W/S` heterozygous, `-`/`N` missing.
  - `bin/mis/plink2genotype.pl` converts PED/MAP into this format (one IUPAC
    code per sample; `0`/`-` are missing, a mix of missing and an allele keeps
    the observed base, the position is the last MAP column). It **streams** the
    PED: memory is about `2 × samples × variants` bytes instead of one Perl
    scalar per field. The default budget is half of `MemAvailable`, which
    normally holds the whole matrix, so one pass is the fastest option and no
    argument is needed; `-maxMem <MiB>` sets a ceiling (the PED is then read in
    blocks, trading time for memory) and `-block <n>` forces the block size.
    The cost of blocking: the PED is re-read once per block, so the overhead
    grows with `variants / block size`; the script warns above 200 passes.
  - Plain or gzip input.
  - With `-SubPop` a header is required: `#CHROM<TAB>POS<TAB>S1 S2 ...` (the
    sample columns are space separated, in any order, matched by name).
  - All-samples path: after the FilterGeno filters (miss/het/biallelic/MAF,
    Cut3base always on), **a missing call becomes a heterozygous
    pseudo-haplotype**; SubPop path: miss/het, strictly biallelic,
    `SeD/(SeD+Max)` MAF, **missing becomes (2,2)**. The two paths differ in
    their filtering and encoding semantics (as they do in the original).
- **`-InVCF` and `-InGenotype` can be given together**: the sites are merged
  into the same per-chromosome table (stable sort, first occurrence wins).
- **`-InPlink` (M7)**: `-InPlink <prefix|*.bed>`, **mutually exclusive** with
  `-InVCF`/`-InGenotype` (PLINK carries its own sample list and genotype
  coding, so the sample order and set cannot be cross-checked; mixing them
  would fail silently).

## PLINK binary input (M7, `-InPlink`)
- Three files with a common prefix: `x.bed` + `x.bim` + `x.fam`; the argument
  takes the prefix `x` or the `x.bed` file itself.
- **`.bed`**: only the standard 2-bit variant-major layout (magic `6c 1b 01`)
  is accepted. Individual-major (`6c 1a`) is refused with a pointer to
  `plink --make-bed`; the mixed SNP+indel mode (`6c 1b 00`, where an indel
  takes two 2-bit hardcalls) is refused with a pointer to dropping the indels.
  The magic is **3 bytes**, and the file size must equal
  `3 + variants × ceil(samples/4)`, otherwise it is an error -- which also
  catches the different layout PLINK uses when a single chromosome has more
  than 524288 variants. Trailing padding bits are ignored when the sample
  count is not a multiple of 4.
  - 2-bit codes: `00` = homozygous A1, `01` = missing, `10` = heterozygous,
    `11` = homozygous A2.
- **`.bim`**: 6 columns `chr / id / cm / bp / a1 / a2`; the position is the
  `bp` column (`cm` is ignored). An allele that is not a single `A/C/G/T`
  (indels, PLINK's `D`/`I`/`0` codes) drops the site through the **Indel
  filter** (like the VCF path dropping a multi-character `REF`/`ALT`);
  `a1 == a2` drops it as non-biallelic.
- **`.fam`**: 6 columns, **column 2 (the IID) is the sample name** (the same
  names `plink --recode vcf` writes into a VCF header); a duplicate IID is an
  error.
- **The filtering and encoding semantics follow `-InVCF`** (not
  `-InGenotype`): the same site filters, the same MAF/Miss/Het denominators
  and the same pseudo-haplotype construction, so **the site set and the
  `NumberPairs` of every distance bin are identical**. What is **not**
  byte-identical to `-InVCF` on real phased data is `r²`/`D'`, because a
  2-bit `.bed` carries no phase:
  - The best/second allele come from the `.bed` counts (`00` counts as the GT
    character `'0'`, `11` as `'1'`, `10` counts once each; on a tie `'0'` is
    the best) -- equivalent to the VCF path scanning the GT characters.
  - The `-Miss`/`-Het` denominator is the **number of samples**; the `-MAF`
    denominator is the **number of observed alleles** (`sed/(sed+max)`).
  - A heterozygote (PLINK has no phase) becomes the pseudo-haplotype
    `(best, sed)`, always encoded as `(0,1)`; a missing call becomes `2`.
  - The known difference against `-InGenotype` (the two old paths differ
    anyway): the genotype path picks best/sed by base letter and divides MAF
    by `2×N` (missing included). Sites whose allele counts are exactly tied and
    sites at the `-MAF` boundary can therefore differ between the two input
    paths (the golden suite records this as an INFO line).
  - **Phase (the key point)**: a `.bed` holds genotypes, not phase, so
    `-InPlink` synthesises a haplotype pair `(best, sed)` per site, while
    `-InVCF` uses the **true phase** in `GT`. Measured on 1000G chr22: the site
    set and the `NumberPairs` of 300,000 bins are identical, but `r²` differs
    by about 0.0055 on average (0.026 at most). Which allele PLINK names `A1`
    does not matter, so this is a phase effect and not an allele-order effect.
    **On unphased data the two are byte-identical**, which is why the older
    golden cases did not expose the difference.
- **`-SubPop`**: the sample names are the `.fam` IIDs (a synthetic
  `#CHROM POS <IID>...` header reuses `resolve_subpop_header` with
  `sample_start=2`); the multi-group and two-column forms work as with VCF.
- **`-L`**: the same as the VCF fallback path -- sites are bucketed by the
  region that owns them and cross-region pairs are not computed.
- **`-OutFilterSNP`**: PLINK has no notion of a "source line", so the passing
  variants are written to `{prefix}.plink.filter.bim` (6-column text in `.bim`
  order, all-samples path only). That mode reads the `.bed` single-threaded, to
  keep the order.
- **`-EHH`**: the same EHH path as the VCF and genotype inputs.
- **Parallelism**: a `.bed` is randomly accessible, so the variant range is
  cut into blocks and decoded in parallel (only when there are at least 4096
  variants, at most 64 blocks). The per-block buffers are concatenated in
  worker order, so the stable sort and the deduplication give the same result
  as the single-threaded run.
- **Known limitation**: the indexed fast path and the position sharding apply
  to `-InVCF` only (see §10/§11 of IMPLEMENTATION.md), so the LD computation
  of a single-chromosome PLINK input runs contig by contig, like
  `-InGenotype` (`-T` still decodes the `.bed` in parallel and still
  distributes several contigs).

## `-L` regions (M3)
- BED: `chr<TAB>start0<TAB>end1` (0-based, half-open `[start,end)`), may be
  `.gz`.
- A single region: `chr` or `chr:start-end` (1-based, inclusive, e.g.
  `chr1:1000-20000`).
- Semantics: only the SNPs inside a region are kept; **LD is computed only
  between SNPs of the same region**; several regions are independent of each
  other and their statistics are aggregated into one output.

## Output: `.stat.gz` (M1/M4, byte-compatible with the original)
```
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
1	0.0020	NA	0.0020	NA	4
2	0.0030	NA	0.0090	NA	3
...
```
- `Dist` = the physical distance in bp, 1..MaxDist(kb)*1000.
- Only bins with `count>0` are written. `Mean=Sum/count`, `%.4f` fixed to 4
  decimals.
- A row whose NumberPairs is 0 does not appear at all.
- With `-OutType 2/3/6/7/8` the D' column carries real values (TF=3 is always
  `0.0000`, because sumD is not accumulated).
- The old `Plot_OnePop.pl / Plot_MultiPop.pl` read it without any change.

## Output: `.stat.gz` (M4, OutType 4/5 histogram)
```
#Dist	R^2	R^2_count
#2SampleSize	120
1	0.00	3
...
```
- OT4 adds the two columns `D'\tD_count` (header
  `#Dist\tR^2\tR^2_count\tD'\tD_count`).
- The `R^2`/`D'` columns hold the bucket representative value `j/100`
  (`setprecision(2)`), `count` is the number of hits, `int(val*100)`; Dist
  starts at 1 and only non-zero buckets are written.
- `#2SampleSize` = 2× the number of effective samples (the subset, with
  SubPop).

## Output: `.vcf.filter.gz` (M3, `-OutFilterSNP`)
- The VCF subset of the sites that pass the filters (Miss/Het/biallelic/MAF),
  in the format of the original.

## Output: `.plink.filter.bim` (M7, `-InPlink -OutFilterSNP`)
- The PLINK variants that pass every filter (Miss/Het/biallelic/MAF), as
  6-column text `chr / id / 0 / bp / a1 / a2`, **in the original `.bim` order,
  without a header**; produced by the all-samples path only (as with the other
  inputs).

## Output: `.genotype.filter.gz` (M6, `-InGenotype -OutFilterSNP`)
- The verbatim subset of the genotype lines that pass the FilterGeno filters
  (Miss/Het/biallelic/MAF + Cut3base); **plain gzip** (the raw bytes equal the
  original ogzstream output, so `cmp` works); the file name is
  `strip_stat_prefix(out)+".genotype.filter.gz"`.
- Produced by the all-samples path only (the SubPop path writes nothing, as in
  the original); no header line.

## Output: `.LD.gz` (M4, OutType 3/6/7/8)
- The pairwise SNP×SNP LD table, one chromosome at a time (plain gzip; the
  file name is `strip_stat_prefix(out)+".LD.gz"`; the row order is chr
  lexicographic → Site1 ascending → Site2 ascending), columns:
  - OT3: `#chr	Site1	Site2	r^2	Dist`
  - OT6: `#chr	Site1	Site2	D'	r^2	Dist`
  - OT7: `#chr	Site1	Site2	D'	LOD	r^2	Dist`
  - OT8: `#chr	Site1	Site2	D'	LOD	r^2	lowCI	highCI	Dist`
- `lowCI/highCI` = the 95% likelihood interval (`low_i/high_i` × 1%, a 101-point
  profile).
- **Every row starts with the chromosome name, and the row length is not
  limited by a buffer**: the original writes these rows through a C++ stream,
  which has no length limit, and so does this implementation (`src/rowfmt.h`)
  -- a row that does not fit is re-formatted into a heap buffer instead of
  being cut off. With a chromosome name of about 90 characters or more the old
  implementation read out of bounds and disagreed with the original; it is
  byte-identical now.

## Output: `.ehh.gz` (M5, `-EHH <chr:Site>`)
- The regional EHH decay around the core site `Site`, with a radius of
  `-MaxDist` (kb, `{prefix}` = `strip_stat_prefix(out)`); **plain gzip** (the
  raw bytes equal the original ogzstream output, so `cmp` works). With `-EHH`
  no `.stat.gz` is written, and the successful path removes
  `{prefix}.stat.gz` again.
- The first two header lines:
  ```
  #Chr	Site	Dist	EHH_all	EHH_0	EHH_1
  #SNP_Number	{count}
  ```
  `count` = the number of sites in the region without a missing call (a site
  with a missing sample is dropped entirely).
- The data rows are the no-missing sites of the region, in ascending pos:
  `chr\tSite\tDis\t{D:.4f}\t{RR:.4f}\t{(D-RR):.4f}`; a site that does not
  reach the `D<0.088` stopping threshold gets a zero row
  `…\t0.0000\t0.0000\t0.0000`.
  - `Dis = Site - core site` (may be negative).
  - `EHH_all = D`, `EHH_0 = RR` (the contribution of the group whose first
    character is `'0'`), `EHH_1 = D-RR`.
- By-products: `{prefix}.ehh.pdf/.ehh.png` (drawn with Rscript, found with
  `which Rscript`, a warning only when it is missing) + `{prefix}.tmp.r`
  (kept with `-OutType 0`, deleted otherwise).
- The data rows start with the chromosome name as well, and their length is
  not limited either (see the `.LD.gz` section).
- A figure that could not be drawn is reported on stderr, but **the exit code
  stays 0** (as in the original).
