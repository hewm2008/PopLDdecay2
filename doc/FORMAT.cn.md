# FORMAT — 输入输出格式

## 输入
- **VCF 4.x**：`-InVCF`，支持 gzip / bgzip / 普通文本；BCF 亦可通过 htslib 透明读取。
- 需 GT 字段（`bcf_get_genotypes`）；无 GT 或相位混杂按旧版语义处理。
- 位点过滤（`-OutFilterSNP` 写出的 `.vcf.filter.gz`）后的文件可直接再喂回 `-InVCF`。
- **InGenotype（M6）**：`-InGenotype`。**原生格式**（每样本一个 IUPAC 记号，与 VCF2Genotype `-WithHeader -NoRef` 输出同风格）：
  ```
  chr1	2310	M A A - A A A A ...
  chr1	2375	R R G R R G - R R ...
  ```
  - `chr<TAB>site<TAB>` 后为每样本一个记号（空格分隔）：`A/C/G/T` 纯合、`M/K/Y/R/W/S` 杂合、`-`/`N` 缺失。
  - `bin/mis/plink2genotype.pl` 输出每样本 **2 个**等位记号，**不兼容**此格式（勿用于 `-InGenotype`）。
  - 支持 gzip / 普通文本输入。
  - `-SubPop` 时需要表头 `#CHROM<TAB>POS<TAB>S1 S2 ...`（样本列空格分隔，顺序任意，按名匹配）。
  - 全样本路径：FilterGeno 过滤（miss/het/双等位/MAF，Cut3base 恒开）后**缺失→杂合伪单倍型**；SubPop 路径：miss/het/恒双等位/`SeD/(SeD+Max)` MAF，**缺失→(2,2)**。两条路径过滤与编码语义不同（与旧版一致）。
- **`-InVCF` 与 `-InGenotype` 可同时给出**：位点合并入同一 per-chr 表（同 chr 稳定排序、去重保首）。

## `-L` 区域（M3）
- BED：`chr<TAB>start0<TAB>end1`（0-based，半开 `[start,end)`），可 `.gz`。
- 单条：`chr` 或 `chr:start-end`（1-based 闭区间，如 `chr1:1000-20000`）。
- 语义：仅保留落在某区间内的 SNP；**LD 仅在属于同一区间的 SNP 对之间计算**；多区间各自独立，统计结果聚合进同一输出。

## 输出：`.stat.gz`（M1/M4，字节兼容旧版）
```
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
1	0.0020	NA	0.0020	NA	4
2	0.0030	NA	0.0090	NA	3
...
```
- `Dist` = 物理距离(bp)，1..MaxDist(kb)*1000。
- 仅输出 `count>0` 的 bin。`Mean=Sum/count`，`%.4f` 固定小数 4 位。
- 若某行 NumberPairs 为 0 则整行不出现。
- `-OutType 2/3/6/7/8` 时 D' 列为实值（TF=3 恒 `0.0000`，因不累加 sumD）。
- 旧 `Plot_OnePop.pl / Plot_MutiPop.pl` 无需修改即可读入。

## 输出：`.stat.gz`（M4，OutType 4/5 直方图）
```
#Dist	R^2	R^2_count
#2SampleSize	120
1	0.00	3
...
```
- OT4 增 `D'\tD_count` 两列（头 `#Dist\tR^2\tR^2_count\tD'\tD_count`）。
- `R^2`/`D'` 列为桶代表值 `j/100`（`setprecision(2)`），`count` 为 `int(val*100)` 命中次数；Dist 从 1 起仅非零桶。
- `#2SampleSize` = 2×有效样本数（SubPop 时按子集）。

## 输出：`.vcf.filter.gz`（M3，`-OutFilterSNP`）
- 通过过滤（Miss/Het/双等位/MAF）的位点的 VCF 子集，格式与旧版一致。

## 输出：`.genotype.filter.gz`（M6，`-InGenotype -OutFilterSNP`）
- 通过 FilterGeno 过滤（Miss/Het/双等位/MAF + Cut3base）的基因型行原文子集；**plain gzip**（原始字节 == 旧版 ogzstream 输出，可 `cmp`）；文件名 = `strip_stat_prefix(out)+".genotype.filter.gz"`。
- 仅全样本路径产生（SubPop 路径不写，同旧版）；无表头行。

## 输出：`.LD.gz`（M4，OutType 3/6/7/8）
- 逐 chr 输出 SNP×SNP 成对 LD 表（plain gzip；文件名 = `strip_stat_prefix(out)+".LD.gz"`；行序 chr 字典序 → Site1 升序 → Site2 升序），列：
  - OT3：`#chr	Site1	Site2	r^2	Dist`
  - OT6：`#chr	Site1	Site2	D'	r^2	Dist`
  - OT7：`#chr	Site1	Site2	D'	LOD	r^2	Dist`
  - OT8：`#chr	Site1	Site2	D'	LOD	r^2	lowCI	highCI	Dist`
- `lowCI/highCI` = 95% 似然区间（`low_i/high_i`×1%，101 点剖面）。

## 输出：`.ehh.gz`（M5，`-EHH <chr:Site>`）
- 以核心位点 `Site` 为中心、半径 = `-MaxDist`（kb，`{prefix}` = `strip_stat_prefix(out)`）的区域 EHH 衰减；**plain gzip**（原始字节 == 旧版 ogzstream 输出，可 `cmp`）。触发 `-EHH` 后不写 `.stat.gz`，成功路径还会删 `{prefix}.stat.gz`。
- 前两行表头：
  ```
  #Chr	Site	Dist	EHH_all	EHH_0	EHH_1
  #SNP_Number	{count}
  ```
  `count` = 区域内无缺失位点数（任一样本缺失的位点整站剔除）。
- 数据行按区域内无缺失位点 pos 升序：`chr\tSite\tDis\t{D:.4f}\t{RR:.4f}\t{(D-RR):.4f}`；未达 `D<0.088` 终止阈值的位点为零值行 `…\t0.0000\t0.0000\t0.0000`。
  - `Dis = Site - 核心Site`（可为负）。
  - `EHH_all = D`、`EHH_0 = RR`（首字符 `'0'` 分组贡献）、`EHH_1 = D-RR`。
- 副产品：`{prefix}.ehh.pdf/.ehh.png`（Rscript 绘图，`which Rscript` 探测，缺失时仅警告）+ `{prefix}.tmp.r`（`-OutType 0` 时保留，否则删除）。