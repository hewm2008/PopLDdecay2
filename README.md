# PopLDdecay2

## Overview

[PopLDdecay2](https://github.com/hewm2008/PopLDdecay2) is a high-performance and memory-efficient successor to
[PopLDdecay](https://github.com/hewm2008/PopLDdecay), designed for large-scale
linkage disequilibrium (LD) decay analysis in modern population genomic datasets.

It preserves the analytical framework, command-line interface, filtering criteria,
and output conventions of the original PopLDdecay while substantially improving
runtime, memory efficiency, and multi-thread scalability. PopLDdecay2 supports
direct analysis of VCF and native genotype formats, as well as efficient
multi-population LD decay analysis within a single run.

On the 1000 Genomes Project Phase 3 chromosome 22 benchmark, PopLDdecay2
achieved a <b>17.5×</b> single-thread speedup with approximately 49% lower peak memory
usage, with further acceleration reaching <b>97.5× at 16 threads and 143× at
32 threads</b> relative to PopLDdecay 3.45.

## Key features

- **High-performance LD decay analysis**  
  Optimized genotype representation, memory organization, and LD computation
  provide substantial acceleration even in single-threaded execution.

- **Scalable multi-thread computation**  
  The `-T` option enables parallel LD calculation across genomic tasks. Indexed
  VCF files additionally support region-based parallel access and position
  sharding for improved scalability.

- **Memory-efficient 2-bit genotype representation**  
  Genotypes are compactly encoded to reduce memory usage and improve cache
  efficiency during pairwise LD calculation.

- **Single-run multi-population analysis**  
  Multiple populations can be analyzed from a single dataset using a two-column
  `sample<TAB>groupid` file, avoiding repeated processing of the same genomic
  input.
  
- **Compatibility with PopLDdecay**  
  PopLDdecay2 preserves the major analytical behavior, filtering criteria,
  command-line conventions, and output formats of the original PopLDdecay.

## Quick start

### 1. Build PopLDdecay2

```bash
./make.sh
```
The executable will be generated at: ```bin/PopLDdecay2```
### 2. Run a basic LD decay analysis using multiple threads
```
bin/PopLDdecay2   -InVCF input.vcf.gz   -OutStat output -T 4
```
This generates: `output.stat.gz`

### 3. Analyze multiple populations in one run
Prepare a two-column sample-to-population file:
```
sample1    AFR
sample2    AFR
sample3    EUR
sample4    EUR
sample5    EAS
```
Then run:
```
bin/PopLDdecay2    -InVCF input.vcf.gz   -SubPop sample.group   -OutStat output   -T 8
```
PopLDdecay2 will generate population-specific LD decay statistics for all groups
defined in `sample.group` within a single run.

### 4. Plot the LD decay curve (same with PopLDdecay1)
For a single populatiom, use:
```
perl bin/Plot_OnePop.pl   -inFile output.stat.gz   -output output
```
For multiple populations, use:
```R
perl bin/Plot_MutiPop.pl    -inList stat.list    -output multi_population
```

<div align="center">  <img src="https://github.com/hewm2008/PopLDdecay2/raw/main/Test/GroupLDecay.png" width="300" alt="GroupLDdecay 示意图"></div>

For all avaiable options, run: `bin/PopLDdecay2 -help`


## Installation
Requires `htslib` (>= 1.15), `zlib`, and `pthreads`.
### Linux
```bash
./make.sh        # auto-detects htslib (>= 1.15); produces bin/PopLDdecay2
# bash  ./make.sh  HTSLIB_ROOT=/path/to/htslib
```
If htslib is installed in a non-standard location, set `HTSLIB_ROOT` before building:
`export HTSLIB_ROOT=/path/to/htslib`
Alternatively, configure CMake directly with:
`cmake -S . -B build -DHTSLIB_ROOT=/path/to/htslib`
### MacOS
```bash
brew install htslib
export HTSLIB_ROOT=$(brew --prefix htslib)
./make.sh        # auto-detects htslib (>= 1.15); produces bin/PopLDdecay2
```

## Main options
| Option | Meaning |
| :--- | :--- |
| `-InVCF` | Input SNP VCF (plain / gzip / bgzip). |
| `-InGenotype` | Input SNP genotype format (see below). |
| `-OutStat` | Output stat file (`{prefix}.stat.gz`). |
| `-SubPop` | Subgroup sample list (repeatable). Single-column file (one sample per line) = one group, labeled by the file basename; two-column `sample<TAB>groupid` = one group per distinct id. One single-column file keeps the legacy `{prefix}.stat.gz` output. |
| `-MaxDist` | Max distance (kb) between two SNPs to consider [300] (capped at 1 Gb − 1 kb). |
| `-MAF` | Min minor allele frequency filter [0.005]. |
| `-Het` | Max het ratio filter [0.88]. |
| `-Miss` | Max missing ratio filter [0.25]. |
| `-EHH` | EHH region decay: `-EHH chr:site` (region = ±`-MaxDist`). |
| `-OutFilterSNP` | Write the passing SNPs (`{prefix}.vcf.filter.gz` / `.genotype.filter.gz`). |
| `-L` | Region restrict: BED / `chr:start-end` / `chr`. |
| `-T` | Threads for parallel LD calculation [4] (capped at 256). With a `.tbi`/`.csi` index and more threads than contigs, large contigs (including a single-chromosome input) are split into overlapping position shards for extra parallelism. |
| `-OutType` | Output type 0-8 (see `-help`). |
| `-help` | Help information |

## Input formats
### VCF format

PopLDdecay2 accepts SNP data in standard VCF format through `-InVCF`.

Supported input types include:

- uncompressed VCF (`.vcf`)
- gzip-compressed VCF (`.vcf.gz`)
- bgzip-compressed and indexed VCF (`.vcf.gz` with `.tbi` or `.csi` index)

Example:

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  S1  S2  S3
chr1    2310 .   A    G    .     PASS    .     GT      0/0 0/1 1/1
```

## Outputs 
PopLDdecay2 generates LD summary statistics by default and can optionally
produce pairwise LD records, filtered variant files, and EHH results depending
on the selected parameters.

| Output file | Description | Generated by |
| --- | --- | --- |
| `{prefix}.stat.gz` | LD decay summary statistics used for downstream plotting and population-level LD decay analysis. | Default output |
| `{prefix}.LD.gz` | Pairwise LD results, including LD statistics for individual SNP pairs. | Selected `-OutType` modes |
| `{prefix}.vcf.filter.gz` | VCF records that pass the specified MAF, missingness, and heterozygosity filters. | `-OutFilterSNP` with `-InVCF` |
| `{prefix}.genotype.filter.gz` | Genotype records that pass the specified filters. | `-OutFilterSNP` with `-InGenotype` |
| `{prefix}.ehh.gz` | Extended haplotype homozygosity (EHH) decay results for the specified genomic site. | `-EHH` |

For multi-population analysis with `-SubPop`, PopLDdecay2 generates
population-specific `.stat.gz` files for each group defined in the
sample-to-population mapping file.

The exact pairwise output format depends on the selected `-OutType` value.
Run:

```bash
bin/PopLDdecay2 -help
## Verification 
bash scripts/golden_diff.sh    # 80 golden cases vs the original binary
```
PopLDdecay2 has been validated against the original PopLDdecay using an 80-case regression test suite covering major input types, analysis options, and output modes.

For identical inputs and parameters:
- decompressed .stat.gz outputs are byte-identical;
- .LD.gz, .vcf.filter.gz, .genotype.filter.gz, and .ehh.gz outputs are byte-identical in compressed form.

All 80 regression cases passed with diff = 0.


## Benchmark
Benchmarks were performed on chr22 of the 1000 Genomes Project Phase 3 dataset (2,504 samples, 1,103,547 SNPs; phased; `-MaxDist 300`) on a shared 80-core RHEL9 x86_64 node.


Thread scaling
| Tool | Threads | Wall time | Peak RSS | vs 3.45 | Byte-identical |
|---|---|---|---|---|---|
| PopLDdecay 3.45 | 1 | 3305 s | 1066 MB | 1.0x | - |
| PopLDdecay2 | 1 | 189 s | 547 MB | **17.5x** | diff = 0 |
| PopLDdecay2 | 16 | 34 s | 712 MB | **97.5x** | diff = 0 |
| PopLDdecay2 | 32 | 23 s | 877 MB | **143x** | diff = 0 |

PopLDdecay2 achieves a 17.5× single-thread speedup while reducing peak memory usage by approximately 49%. Additional parallelization further reduces runtime to 34 s with 16 threads and 23 s with 32 threads.

Decompressed .stat.gz outputs were byte-identical to those generated by PopLDdecay 3.45 at all tested thread counts.

**Multi-population analysis**
PopLDdecay2 can analyze multiple populations in a single run using a two-column `sample<TAB>groupid` file.
| Method | Populations | Threads | Wall time | Peak RSS | Speedup |
|---|---:|---:|---:|---:|---:|
| PopLDdecay 3.45 | 5 sequential runs | 1 | 2720.79 s | — | 1.0× |
| PopLDdecay2 | 5 populations in one run | 8 | 57.57 s | 616 MB | **47.3×** |

The benchmark included the five 1000 Genomes super-populations: AFR, AMR, EAS, EUR, and SAS. Population-specific decompressed .stat.gz outputs were byte-identical to the corresponding PopLDdecay 3.45 results.
**Reproduce the benchmarks**

Benchmark commands: `bash scripts/bench_chr22.sh`, `bash scripts/bench_subpop.sh`

## Citation

If you use **PopLDdecay2**, please cite:
> **[[PopLDdecay2](https://github.com/hewm2008/PopLDdecay2)](https://github.com/hewm2008/PopLDdecay2)**
## License

PopLDdecay2 is open-source software released under the **MIT License**.  
You are free to use, modify, and distribute the software in accordance with the terms of the license.

See the [LICENSE](LICENSE) file for details.
## Discussion
------------
- [:email:](https://github.com/hewm2008/PopLDdecay2) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>

######################swimming in the sky and flying in the sea #############################
