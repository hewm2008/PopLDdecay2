# PopLDdecay2

<b>PopLDdecay2：超快、省内存的连锁不平衡（LD）衰减分析工具</b>

[English](README.md) | [输入输出格式](doc/FORMAT.cn.md) | [复现脚本](RunTest/README.txt)

**版本：v0.80**

PopLDdecay2 算出的结果与 [PopLDdecay](https://github.com/hewm2008/PopLDdecay) (hewm2008) 完全相同——命令行相同、过滤标准相同、**输出逐字节一致**——但快得多、内存约省一半。在 1000 Genomes Phase 3 的 22 号染色体上，单线程加速 **17.6×**，32 线程 **147×**。

PopLDdecay2：可以先引这儿 [hewm2008/PopLDdecay2](https://github.com/hewm2008/PopLDdecay2)，为了让更多人知晓,辛苦尽量引用。

### 1) 安装与参数
------------

* 1) 安装

`sh make.sh` 即可编译，可执行文件在 <b>`bin/PopLDdecay2`</b>。
  </br> <b>Linux / Unix</b> 与 <b>macOS</b>
  <pre>
        tar -zxvf  PopLDdecay2XXX.tar.gz
        cd PopLDdecay2XXX;
        sh make.sh;                             # 需要 htslib (>= 1.18)、zlib、pthreads
        ./bin/PopLDdecay2 -help
  </pre>

**环境要求：** 支持 C++17 的编译器（GCC 8+ 或同类编译器）

htslib 不在标准位置时：`sh make.sh HTSLIB_ROOT=/path/to/htslib`，或 `cmake -S . -B build -DHTSLIB_ROOT=/path/to/htslib`。
  </br> macOS：`brew install htslib` 并 `export HTSLIB_ROOT=$(brew --prefix htslib)`。

编译在 x86 上会显式请求硬件 POPCNT 指令（内核用 popcount 统计单倍型），并关闭 FP 收缩，使输出不依赖编译器的默认 target。两者都可关闭：`sh make.sh ENABLE_MPOPCNT=0 ENABLE_FP_CONTRACT_OFF=0`。

**编译开关**（都以 `sh make.sh KEY=VALUE` 的形式给出）：

| 开关 | 作用 |
| --- | --- |
| `HTSLIB_ROOT=<prefix>` | htslib 前缀；需含 `include/htslib` 与 `lib/libhts` |
| `CXX=<compiler>` | 用这个编译器编译，不使用自动探测到的那一个 |
| `NATIVE=1` | `-march=native`：更快，但与机器绑定；边界处的 1-ULP 舍入可能让与黄金摘要的逐字节比对失败——要比对时请保持 `0` |
| `LTO=1` | 链接期优化：语义不变，输出逐字节一致 |
| `KEEP_BUILD=1` | 保留中间产物 `build/`，不删除 |
| `ENABLE_MPOPCNT=0` | 不请求硬件 POPCNT 指令 |
| `ENABLE_FP_CONTRACT_OFF=0` | 不关闭浮点收缩 |

后两个只影响速度，不影响输出字节。

* 2) 参数

```sh
	Usage: PopLDdecay2 -InVCF  <in.vcf.gz>  -OutStat <out.stat>

		-InVCF         <str>     Input SNP VCF Format
		-InGenotype    <str>     Input SNP Genotype Format
		-InPlink       <str>     Input PLINK [bed+bim+fam] file prefix
		-OutStat       <str>     OutPut Stat Dist ~ r^2/D' File

		-SubPop        <str>     SubGroup Sample File List (repeatable)
		-MaxDist       <int>     Max Distance (kb) between two SNP [300]
		-MAF           <float>   Min minor allele frequency filter [0.005]
		-Het           <float>   Max ratio of het allele filter [0.88]
		-Miss          <float>   Max ratio of miss allele filter [0.25]
		-EHH           <str>     To Run EHH Region decay set StartSite [NA]
		-OutFilterSNP            OutPut the final SNP to calculate
		-L             <str>     Region restrict: BED/chr:start-end/chr
		-T             <int>     Threads for parallel LD calculation [4]
		-OutType       <int>     OutType [0-8], see -help
		-help                    Show more help [hewm2008]
```

更多说明请用 <b>-help</b>，并看 [doc/FORMAT.cn.md](doc/FORMAT.cn.md)。

数值按严格方式解析：`-MAF`/`-Het`/`-Miss` 必须是有限小数且落在 `[0, 1]`，`-MaxDist`/`-OutType`/`-T` 必须是整数，`-EHH` 的坐标必须是整数。写错时会带着参数名与原值报错，而不是被静默改写（过去 `-MAF abc` 会变成 0，即 MAF 过滤失效）。边界值 `-MAF 0 -Miss 1 -Het 1` 仍然合法。

* 3) 输入

| 格式 | 扩展名 | 索引快路径 |
| :--- | :--- | :--- |
| VCF | `.vcf` | 无 |
| VCF，gzip | `.vcf.gz`（普通 gzip） | 无 |
| VCF，bgzip | `.vcf.gz`（BGZF） | 有，配 `.tbi`/`.csi` |
| 基因型 | `-InGenotype`，`.gz` 或普通 | 不适用 |
| PLINK 二进制 | `-InPlink prefix(bed/bim/fam)`| 不适用 |

普通 **gzip** 压缩的 VCF 读取同样很快，但无法建索引，因此 `-T > 1` 会退回顺序读（LD 计算仍并行）。要启用索引并行读，先 bgzip 再建索引：`bgzip -@ 8 in.vcf && tabix -@ 8 in.vcf.gz`。

**基因型**（`-InGenotype`）：`chr<TAB>site<TAB>` 后每个样本一个 IUPAC 记号（`A/C/G/T` 纯合、`M/K/Y/R/W/S` 杂合、`-`/`N` 缺失）。

**PLINK**（`-InPlink`）：直接读 `.bed`/`.bim`/`.fam` 三件套，无需转换。样本名取 `.fam` 的 ID，与 `plink --recode vcf` 写进 VCF 头的名字一致，`-SubPop` 名单可原样迁移。

三条路径并不完全相同，比较结果时需注意：

| | `-InVCF` | `-InGenotype` | `-InPlink` |
| :--- | :--- | :--- | :--- |
| 相位 | 用 `GT` 里的真实相位 | 没有 | **`.bed` 不存相位** |
| `-MAF` 分母 | 观测到的等位数 | `2 x 样本数` | 观测到的等位数 |
| 与 3.45 逐字节一致 | 是 | 是 | 3.45 没有 PLINK 输入 |

`-InVCF` 与 `-InPlink` 选中的**位点相同**、**配对数相同**，但真实定相数据上 `r^2` 略有差异（1000 Genomes chr22 上平均约 0.005），因为 2-bit 的 `.bed` 装不下 `GT` 里的相位。

* 4) 输出

`-OutStat` 是**前缀**而不是文件名：统计结果写到 `<prefix>.stat.gz`。值末尾的 `.stat.gz` 或 `.stat` 会先被剥掉，所以 `-OutStat chr22`、`-OutStat chr22.stat` 与 `-OutStat chr22.stat.gz` 最终都得到 `chr22.stat.gz`。带 `-SubPop` 时每个群体各写一个 `<prefix>.<label>.stat.gz`。

| 模式 | 写出的文件 |
| --- | --- |
| 默认（`-OutType 1`） | `<prefix>.stat.gz`：按距离分箱的 r^2 均值 |
| `-OutType 2` | `<prefix>.stat.gz`：r^2 与 D' 的汇总 |
| `-OutType 3`、`6`、`7`、`8` | 上述统计，外加逐对 LD 值 `<prefix>.LD.gz` |
| `-OutType 4`、`5` | 供分布图使用的直方图式统计 |
| `-EHH chr:site` | `<prefix>.ehh.gz`；若能找到 `Rscript`，另出 `<prefix>.ehh.pdf` / `<prefix>.ehh.png` |
| `-OutFilterSNP` | 通过过滤的位点，格式随输入而定：`<prefix>.vcf.filter.gz`、`<prefix>.genotype.filter.gz` 或 `<prefix>.plink.filter.bim` |

**覆盖语义。** 目标 `<prefix>.stat.gz` 会在启动时被以写方式打开，用来检查是否可写——这会**截断同名的旧结果**，哪怕这次运行随后因输入出错退出。每个群体的 `<prefix>.<label>.stat.gz` 同样会被截断。PopLDdecay 3.45 的做法完全一样（它也是用 `ogzstream` 打开同一个文件做同样的检查），所以这是有意复刻的行为，不是疏漏：需要保留旧结果时请换一个前缀。其余文件（`.LD.gz`、`.ehh.gz`、过滤文件）只在对应功能启用时才写，且写失败不会留下半截文件。

* 5) 常见问题

| 现象 | 处理办法 |
| --- | --- |
| 编译时报 "htslib not found" | 装好 htslib，并把 `HTSLIB_ROOT` 指到它的前缀；该前缀需含 `include/htslib` 与 `lib/libhts` |
| 编译时报 "htslib 1.16 or newer is required" | 升级 htslib；1.18 或更新更快 |
| `-T > 1` 仍然串行读 VCF | 普通 gzip 的 VCF 无法建索引：改用 `bgzip -@ 8 in.vcf && tabix -@ 8 in.vcf.gz` |
| `-SubPop` 里的某个群体找不到 | 样本名必须与输入表头完全一致：VCF 的 `#CHROM` 行、基因型文件的表头、PLINK `.fam` 的 ID |
| PLINK 路径与 VCF 路径的 r^2 对不上 | `.bed` 存的是基因型而不是相位，相位数据上两条路径本来就会有差异（见上表） |
| 需要每一对位点的 LD 值 | 用 `-OutType 3`、`6`、`7` 或 `8`；各列定义见 [doc/FORMAT.cn.md](doc/FORMAT.cn.md) |
| 某个参数被拒绝并报错 | 数值按严格方式解析：`-MAF`/`-Het`/`-Miss` 必须是 `[0, 1]` 内的有限小数，`-MaxDist`/`-OutType`/`-T` 必须是整数（见上文） |

### 2) 示例
------------

包内自带一份小数据集（50,000 个双等位 SNP × 2,504 样本，19 MB），因此整条流程不联网也能走通。

* 1) 单个群体的 LD 衰减

```sh
# 1.1) 运行
        ./bin/PopLDdecay2  -InVCF RunTest/00.GetData/base_chr22_50k.vcf.gz \
                           -OutStat chr22 -MaxDist 300 -T 8

# 1.2) 看结果：距离分箱
        zcat chr22.stat.gz | head -5
        # #Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
        # 6	0.4153	NA	2.9073	NA	7
        # 7	0.8607	NA	3.4428	NA	4

# 1.3) 画图（需要 R）
        perl bin/Plot_OnePop.pl  -inFile chr22.stat.gz  -output chr22
```

  </br>`Dist` 是碱基（bp）；没有配对的距离不会出现在文件里；`NA` 表示该模式不填这一列。

* 2) 一次运行分析五个群体

```sh
# 2.1) 两列文件：样本 <TAB> 群体（已提供：RunTest/04.Fig3Run/sample.group）

# 2.2) 一次运行，每个群体一个 .stat.gz
        ./bin/PopLDdecay2  -InVCF RunTest/00.GetData/base_chr22_50k.vcf.gz \
                           -SubPop RunTest/04.Fig3Run/sample.group -OutStat mp -T 8
        # -> mp.AFR.stat.gz、mp.AMR.stat.gz、mp.EAS.stat.gz、mp.EUR.stat.gz、mp.SAS.stat.gz

# 2.3) 画到同一张图
        printf 'mp.AFR.stat.gz\tAFR\nmp.AMR.stat.gz\tAMR\nmp.EAS.stat.gz\tEAS\nmp.EUR.stat.gz\tEUR\nmp.SAS.stat.gz\tSAS\n' > pops.list
        perl bin/Plot_MultiPop.pl  -inList pops.list  -output mp
```

  </br>`-SubPop` 的样本名在输入文件自己的表头里匹配（VCF 的 `#CHROM` 行、基因型文件的 `#CHROM POS S1 S2 ...` 表头、PLINK `.fam` 的 ID）。名字找不到会提示，重复会报错，群体样本少于 3 个会报错。`-L` 可给 BED 文件（0-based、半开）或 `chr:start-end`（1-based、闭区间）；**LD 只在属于同一区域的 SNP 之间计算**。

* 3) 运行过程中的提示

```sh
[PopLDdecay2] input VCF /path/base_chr22_50k.vcf.gz (indexed)
[PopLDdecay2] filters MAF=0.005 Het=0.88 Miss=0.25 MaxDist=300kb
[PopLDdecay2] threads=8 groups=1 (all samples)
[PopLDdecay2] output chr22.stat.gz
[PopLDdecay2] sites read=52469 kept=52469 (indel=0 multiallelic=0 filters=0)
[PopLDdecay2] wrote chr22.stat.gz
```

  </br>这些信息走 **stderr**；stdout 保持与旧版完全一致：

```sh
Used [perl  ../bin/Plot_XX.pl ] to Plot the LDdecay
```

  </br>长任务在 LD 计算阶段每五分钟追加一行进度，写明染色体、已统计的配对数与已用时间：

```sh
[PopLDdecay2] 22 pairs=19494095 elapsed=300s
```

### 3) 性能基准
------------

1000 Genomes Phase 3 的 22 号染色体（2,504 样本、1,103,547 条记录、已定相、bgzip + `.tbi`），`-MaxDist 300 -MAF 0.005`，共享的 80 核 RHEL9 x86_64 节点：

| 工具 | 线程数 | 墙钟时间 | 峰值 RSS | vs 3.45 | 字节一致 |
|---|---|---|---|---|---|
| PopLDdecay 3.45 | 1 | 3339 s | 1066 MB | 1.0× | - |
| PopLDdecay2 | 1 | 189 s | 547 MB | **17.6×** | diff = 0 |
| PopLDdecay2 | 16 | 33 s | 714 MB | **102×** | diff = 0 |
| PopLDdecay2 | 32 | 23 s | 878 MB | **147×** | diff = 0 |

  </br>（以上是论文中的数字，测量时节点上有其它作业。v0.80 的内核比测量时用的那一版快约一倍。）

1103547 条记录里只有 217103 个真正进入 LD 计算：44030 个 indel 与 7835 个多位点被跳过，其余又被 `-MAF 0.005` 去掉约 80%。

**主要优化：**
- <b>2-bit 基因型</b> —— 单倍型按位打包，一对位点的 2×2 表用 popcount 统计，而不是逐基因型拼装。
- <b>每词一次 popcount</b> —— 一对位点都没有缺失时（1000 Genomes 里 96% 的位点如此），2×2 表由每词一次 popcount 加上按位点的等位计数得出。
- <b>索引并行读</b> —— 带 tabix 索引的 bgzip VCF 按 contig 分配 worker；单条染色体会被拆成重叠的位置分片，于是"读一条 contig"与"算另一条"重叠。
- <b>硬件 popcount</b> —— 编译时请求 POPCNT 指令，比严格 x86-64 target 生成的软件实现快 3.4 倍。

`-T` 带来两样东西，都需要有活可干：并行读需要 tabix 索引，摊开配对需要线程数多于 contig 数。内存随 `-MaxDist`（距离分箱）与 `-T`（每 worker 一份分箱）增长；成对模式（`-OutType 3/6/7/8`）改为吃磁盘。

### 4) 简介
------------

LD 衰减分析测量的是两个位点之间的等位相关性随物理距离的下降。对窗口内的每一对 SNP 计算 r²（以及归一化的 D' 系数），再按距离分箱：

```sh
                r2(i,j) = (P(Ai*Aj) - P(Ai)*P(Aj)) ^ 2 / ( P(Ai)*(1-P(Ai)) * P(Aj)*(1-P(Aj)) )
```

  </br>曲线的形状反映这个群体经历了多少重组，也是群体结构推断、基于 LD 的亲缘估计以及关联研究设计的基础。下面是结果长什么样：1000 Genomes 的五个超级群体，一次遍历一个文件（`RunTest/04.Fig3Run`）：

![example1.png](doc/figures/group_lddecay.png)

### 5) 结果
------------

成对数值不只是自洽，而且是对的。`RunTest/08.Accuracy` 取 50k 切片里 MAF ≥ 0.1 的前 100 个位点（2,504 样本、已定相），把 300 kb 窗口内的 3,590 对位点与"直接数 VCF 里的定相单倍型"这个真值、以及 PLINK 2.0 和 Haploview 4.2 作比较：

| | R^2 对真值 | Pearson | mean\|d\| | max\|d\| |
|---|---|---|---|---|
| **PopLDdecay2** | | **1.000000** | 0.000025 | 0.000050 |
| Haploview 4.2 | | 0.999999 | 0.000238 | 0.000500 |
| PLINK 2.0 | | 0.979098 | 0.020092 | 0.521827 |

  </br>PopLDdecay2 与 Haploview 给出的 r² 相同（只差打印位数），且都与真值一致。PLINK 1.9 与 2.0 彼此一致，却与真值不符：它们**不使用 VCF 里已有的相位**，而是从二倍体计数重新估计单倍型频率。`RunTest/08.Accuracy/step4a` 直接证实了这点——同一批位点，PLINK 对定相版与未定相版给出完全相同的数，而 PopLDdecay2 用到了相位。

相同输入与参数下，所有数据输出与原版 PopLDdecay 逐字节一致：

```sh
OLD=/path/to/PopLDdecay-3.45/bin/PopLDdecay bash RunTest/06.Verify/golden_diff.sh
```

  </br>CI 门禁既不需要那个二进制也不需要下载：`bash RunTest/06.Verify/step1_checkOutput.sh`（108 个用例：88 个输出摘要、17 个行为断言、3 个文档检查）。

### 6) 交流
------------
- [:email:](https://github.com/hewm2008/PopLDdecay) hewm2008@gmail.com / hewm2008@qq.com
- QQ 群：125293663

