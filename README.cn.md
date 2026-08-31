# PopLDdecay2

## 概览

[PopLDdecay2](https://github.com/hewm2008/PopLDdecay2)是 [PopLDdecay](https://github.com/hewm2008/PopLDdecay) 的高性能、内存高效的继任版本，专为现代群体基因组规模数据的连锁不平衡（LD）衰减分析而设计。

它保留原版 PopLDdecay 的分析框架、命令行接口、过滤标准和输出约定，同时显著提升了运行时间、内存效率与多线程可扩展性。PopLDdecay2 支持直接分析 VCF 与原生基因型格式，并可在单次运行内高效进行多群体 LD 衰减分析。

在 1000 Genomes Project Phase 3 的 22 号染色体基准测试中，PopLDdecay2 单线程加速达 **17.5×**，峰值内存占用降低约 49%；进一步并行化后，相对 PopLDdecay 3.45 分别在 16 线程达到 **97.5×**、32 线程达到 **143×** 的加速。

## 核心特性

- **高性能 LD 衰减分析**  
  优化后的基因型表示、内存组织与 LD 计算即使在单线程执行下也能提供大幅加速。

- **可扩展的多线程计算**  
  `-T` 选项可在基因组任务之间并行计算 LD。带索引的 VCF 文件额外支持基于区域的并行访问以及位置分片，进一步提升可扩展性。

- **内存高效的 2-bit 基因型表示**  
  基因型被紧凑编码，以减少内存占用并提升成对 LD 计算时的缓存效率。

- **单次运行多群体分析**  
  通过两列 `sample<TAB>groupid` 文件，可在同一数据集上分析多个群体，避免重复处理同一份基因组输入。

- **与 PopLDdecay 兼容**  
  PopLDdecay2 保留了原版 PopLDdecay 的主要分析行为、过滤标准、命令行约定与输出格式。

## 快速开始

### 1. 构建 PopLDdecay2

```bash
./make.sh
```
可执行文件将生成于：`bin/PopLDdecay2`
### 2. 使用多线程运行基础 LD 衰减分析
```
bin/PopLDdecay2   -InVCF input.vcf.gz   -OutStat output -T 4
```
这将生成：`output.stat.gz`

### 3. 单次运行分析多个群体
准备两列式的样本到群体映射文件：
```
sample1    AFR
sample2    AFR
sample3    EUR
sample4    EUR
sample5    EAS
```
然后运行：
```
bin/PopLDdecay2    -InVCF input.vcf.gz   -SubPop sample.group   -OutStat output   -T 8
```
PopLDdecay2 将在单次运行内为 `sample.group` 中定义的所有群体生成群体特异的 LD 衰减统计结果。

### 4. 绘制 LD 衰减曲线（与 PopLDdecay1 相同）
单群体使用：
```
perl bin/Plot_OnePop.pl   -inFile output.stat.gz   -output output
```
多群体使用：
```
perl bin/Plot_MutiPop.pl    -inList stat.list    -output multi_population
```
<div align="center">  <img src="https://github.com/hewm2008/PopLDdecay2/raw/main/Test/GroupLDecay.png" width="300" alt="GroupLDdecay 示意图"></div>
查看所有可用选项，运行：`bin/PopLDdecay2 -help`


## 安装
需要 `htslib`（>= 1.15）、`zlib` 与 `pthreads`。
### Linux
```bash
./make.sh        # 自动检测 htslib（>= 1.15）；生成 bin/PopLDdecay2
# bash  ./make.sh  HTSLIB_ROOT=/path/to/htslib
```
如果 htslib 安装在非标准位置，请在编译前设置 `HTSLIB_ROOT`：
`export HTSLIB_ROOT=/path/to/htslib`
或者直接使用 CMake 配置：
`cmake -S . -B build -DHTSLIB_ROOT=/path/to/htslib`
### MacOS
```bash
brew install htslib
export HTSLIB_ROOT=$(brew --prefix htslib)
./make.sh        # 自动检测 htslib（>= 1.15）；生成 bin/PopLDdecay2
```

## 主要选项
| 选项 | 含义 |
| :--- | :--- |
| `-InVCF` | 输入 SNP VCF（普通 / gzip / bgzip）。 |
| `-InGenotype` | 输入 SNP 基因型格式（见下）。 |
| `-OutStat` | 输出统计文件（`{prefix}.stat.gz`）。 |
| `-SubPop` | 子群体样本列表（可重复）。单列文件（每行一个样本）= 一个群体，以文件名 basename 标记；两列 `sample<TAB>groupid` = 每个不同 id 一个群体。单个单列文件保持旧版 `{prefix}.stat.gz` 输出。 |
| `-MaxDist` | 两个 SNP 之间的最大距离（kb）[300]（上限为 1 Gb − 1 kb）。 |
| `-MAF` | 最小次要等位基因频率过滤 [0.005]。 |
| `-Het` | 最大杂合度比例过滤 [0.88]。 |
| `-Miss` | 最大缺失比例过滤 [0.25]。 |
| `-EHH` | EHH 区域衰减：`-EHH chr:site`（区域 = ±`-MaxDist`）。 |
| `-OutFilterSNP` | 写入通过的 SNP（`{prefix}.vcf.filter.gz` / `.genotype.filter.gz`）。 |
| `-L` | 区域限制：BED / `chr:start-end` / `chr`。 |
| `-T` | 用于并行 LD 计算的线程数 [4]（上限 256）。配合 `.tbi`/`.csi` 索引，当线程数多于 contig 数时，大 contig（包括单染色体输入）会被拆分为重叠的位置分片以获取额外并行度。 |
| `-OutType` | 输出类型 0-8（见 `-help`）。 |
| `-help` | 帮助信息 |

## 输入格式
### VCF 格式

PopLDdecay2 通过 `-InVCF` 接受标准 VCF 格式的 SNP 数据。

支持的输入类型包括：

- 未压缩的 VCF（`.vcf`）
- gzip 压缩的 VCF（`.vcf.gz`）
- bgzip 压缩并带索引的 VCF（`.vcf.gz` 并带 `.tbi` 或 `.csi` 索引）

示例：

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  S1  S2  S3
chr1    2310 .   A    G    .     PASS    .     GT      0/0 0/1 1/1
```

## 输出
PopLDdecay2 默认生成 LD 汇总统计，并可根据所选参数，按需输出成对 LD 记录、过滤后的变异文件以及 EHH 结果。

| 输出文件 | 描述 | 生成方式 |
| --- | --- | --- |
| `{prefix}.stat.gz` | LD 衰减汇总统计，用于下游绘图与群体级 LD 衰减分析。 | 默认输出 |
| `{prefix}.LD.gz` | 成对 LD 结果，包含单个 SNP 对的 LD 统计。 | `-OutType` 选择的模式 |
| `{prefix}.vcf.filter.gz` | 通过指定 MAF、缺失率与杂合度过滤的 VCF 记录。 | `-OutFilterSNP` 搭配 `-InVCF` |
| `{prefix}.genotype.filter.gz` | 通过指定过滤的基因型记录。 | `-OutFilterSNP` 搭配 `-InGenotype` |
| `{prefix}.ehh.gz` | 指定基因组位点的扩展单倍型纯合子（EHH）衰减结果。 | `-EHH` |

对于使用 `-SubPop` 的多群体分析，PopLDdecay2 为样本到群体映射文件中定义的每个群体生成群体特异的 `.stat.gz` 文件。

成对输出的具体格式取决于所选的 `-OutType` 值。运行：

```bash
bin/PopLDdecay2 -help
## 校验
bash scripts/golden_diff.sh    # 80 个 golden 用例对比原始二进制
```
PopLDdecay2 已使用覆盖主要输入类型、分析选项与输出模式的 80 个回归用例，与原版 PopLDdecay 进行了验证。

对于相同的输入与参数：
- 解压后的 .stat.gz 输出字节级一致；
- .LD.gz、.vcf.filter.gz、.genotype.filter.gz 与 .ehh.gz 输出在压缩形式上字节级一致。

全部 80 个回归用例均通过，diff = 0。


## 基准测试
基准测试在 1000 Genomes Project Phase 3 数据集的 chr22 上进行（2504 个样本，1103547 个 SNP；已定相；`-MaxDist 300`），运行于一台共享的 80 核 RHEL9 x86_64 节点。


线程扩展
| 工具 | 线程数 | 墙钟时间 | 峰值 RSS | vs 3.45 | 字节一致 |
|---|---|---|---|---|---|
| PopLDdecay 3.45 | 1 | 3305 s | 1066 MB | 1.0x | - |
| PopLDdecay2 | 1 | 189 s | 547 MB | **17.5x** | diff = 0 |
| PopLDdecay2 | 16 | 34 s | 712 MB | **97.5x** | diff = 0 |
| PopLDdecay2 | 32 | 23 s | 877 MB | **143x** | diff = 0 |

PopLDdecay2 单线程加速达 17.5×，同时峰值内存占用降低约 49%。进一步并行化后，运行时间分别降至 16 线程的 34 s 与 32 线程的 23 s。

在所有测试线程数下，解压后的 .stat.gz 输出均与 PopLDdecay 3.45 生成的输出字节级一致。

**多群体分析**
PopLDdecay2 可使用两列 `sample<TAB>groupid` 文件单次运行分析多个群体。
| 方法 | 群体数 | 线程 | 墙钟时间 | 峰值 RSS | 加速 |
|---|---|---:|---:|---:|---:|
| PopLDdecay 3.45 | 5 次顺序运行 | 1 | 2720.79 s | — | 1.0× |
| PopLDdecay2 | 单次运行 5 个群体 | 8 | 57.57 s | 616 MB | **47.3×** |

该基准测试包含五个 1000 Genomes 超级群体：AFR、AMR、EAS、EUR 与 SAS。群体特异的解压 .stat.gz 输出与相应的 PopLDdecay 3.45 结果字节级一致。
**复现基准测试**

基准测试命令：`bash scripts/bench_chr22.sh`、`bash scripts/bench_subpop.sh`

## 引用

如果您使用 **PopLDdecay2**，请引用：
> **[[PopLDdecay2](https://github.com/hewm2008/PopLDdecay2)](https://github.com/hewm2008/PopLDdecay2)**
## 许可证

PopLDdecay2 是根据 **MIT License** 发布的开源软件。
您可以根据许可证条款自由使用、修改和分发本软件。

详情请参见 [LICENSE](LICENSE) 文件。
## 交流
------------
- [:email:](https://github.com/hewm2008/PopLDdecay2) hewm2008@gmail.com / hewm2008@qq.com
- 加入<b><i> QQ 群：125293663</b></i>

######################swimming in the sky and flying in the sea #############################
