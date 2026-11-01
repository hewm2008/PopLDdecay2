# FORMAT — 输入输出格式

## 命令行参数的取值范围
| 参数 | 取值 | 越界/非法时的行为 |
|---|---|---|
| `-MAF` / `-Het` / `-Miss` | 有限小数，`[0, 1]`（闭区间） | 报错退出（退出码 1），信息含参数名与原值 |
| `-MaxDist` | 整数（kb） | 非整数报错退出；`< 1` 仍按旧行为重置为 1 并提示；过大（> 1 Gb−1 kb）截断并提示 |
| `-OutType` | 整数 `0`–`8` | 越界或非整数报错退出 |
| `-Method` | 整数，当前仅 `1` | 其他值报错退出（文案同旧版：`Method 2 not implemented yet`） |
| `-T` / `-Threads` | 整数，`1`–`256` | 非整数报错退出；`< 1` 重置为 1、`> 256` 截断为 256，均提示 |
| `-EHH` | `chr:整数` | 缺 `:` 或坐标非整数报错退出 |

- **数值必须是完整数字**：`atoi`/`atof` 无法区分 `"abc"` 与 `0`，旧版据此把 `-MAF abc` 当作 0（MAF 过滤失效）、`-OutType abc` 当作 0（只输出表头）、`-MaxDist 1e9` 当作 1 kb、`-EHH chr1:abc` 把核心位点当 0。这些现在一律报错。
- 边界值仍然合法：`-MAF 0 -Miss 1 -Het 1`（三个过滤全放开）、`-MaxDist 0`（重置为 1）。
- **缺参判定**：若下一个记号是另一个已知选项名，视为缺参（`-OutStat -MaxDist 300` 不再把输出文件命名为 `-MaxDist`）；仅以 `-` 开头但不是选项名的值仍按值接受（路径可以是 `-foo`）。
- 与旧版的关系：旧版同样用 `atoi`/`atof`，因此也会接受这些写法。差异**只在错误路径**（报错信息与退出码），成功运行的输出字节不变。`-EHH` 失败时的退出码仍为 0（与旧版一致），该约束由 CI 用例钉住。

## 输出文件的覆盖语义
- 启动时会以写入方式打开目标 `{prefix}.stat.gz` 以检查可写性——**这会截断同名旧结果**，即使随后因输入错误而失败。旧版 3.45 行为完全相同（它用 `ogzstream` 打开同一文件做同样的检查），此处是**刻意复刻**，不是疏漏。
- 同理，多群体（`-SubPop`）时每个群体的 `{prefix}.{label}.stat.gz` 都会被截断。
- `.LD.gz` / `.ehh.gz` / `.vcf.filter.gz` / `.genotype.filter.gz` / `.plink.filter.bim` 只在对应功能开启时写出；写出失败时不会留下残缺文件（失败即删除并返回非 0 退出码）。

## 输入
- **VCF**：`-InVCF`，支持 gzip / bgzip / 普通文本。**不支持 BCF**：读取器按文本行读取文件。
- 需 GT 字段（`bcf_get_genotypes`）；无 GT 或相位混杂按旧版语义处理。
- 位点过滤（`-OutFilterSNP` 写出的 `.vcf.filter.gz`）后的文件可直接再喂回 `-InVCF`。
- **InGenotype（M6）**：`-InGenotype`。**原生格式**（每样本一个 IUPAC 记号，与 VCF2Genotype `-WithHeader -NoRef` 输出同风格）：
  ```
  chr1	2310	M A A - A A A A ...
  chr1	2375	R R G R R G - R R ...
  ```
  - `chr<TAB>site<TAB>` 后为每样本一个记号（空格分隔）：`A/C/G/T` 纯合、`M/K/Y/R/W/S` 杂合、`-`/`N` 缺失。
  - `bin/mis/plink2genotype.pl` 把 PED/MAP 转为此格式（每样本 1 个 IUPAC 记号；`0`/`-` 记为缺失，缺失+等位混合保留已观测碱基，位置取 MAP 最后一列）。**流式处理**：内存约为 `2 × 样本数 × 变异数` 字节，不再把整个 PED 按字段展开为 Perl 标量。默认预算取 `MemAvailable` 的一半，通常一趟读完即最优，无需传参；`-maxMem <MiB>` 设内存上限（超限则按块处理，以时间换内存），`-block <n>` 可强制块大小。
    分块代价：PED 每块重读一遍，总开销随 `变异数 / 块大小` 增长；趟数超过 200 时脚本会告警。
  - 支持 gzip / 普通文本输入。
  - `-SubPop` 时需要表头 `#CHROM<TAB>POS<TAB>S1 S2 ...`（样本列空格分隔，顺序任意，按名匹配）。
  - 全样本路径：FilterGeno 过滤（miss/het/双等位/MAF，Cut3base 恒开）后**缺失→杂合伪单倍型**；SubPop 路径：miss/het/恒双等位/`SeD/(SeD+Max)` MAF，**缺失→(2,2)**。两条路径过滤与编码语义不同（与旧版一致）。
- **`-InVCF` 与 `-InGenotype` 可同时给出**：位点合并入同一 per-chr 表（同 chr 稳定排序、去重保首）。
- **`-InPlink`（M7）**：`-InPlink <prefix|*.bed>`，与 `-InVCF`/`-InGenotype` **互斥**（PLINK 自带样本表与基因型编码，样本顺序/集合无法交叉校验，混用会静默出错）。

## PLINK 二进制输入（M7，`-InPlink`）
- 三个文件同前缀：`x.bed` + `x.bim` + `x.fam`；参数可给前缀 `x` 或直接给 `x.bed`。
- **`.bed`**：只接受 variant-major（magic `6c 1b 01`）**标准 2-bit 布局**；个体优先（`6c 1a`）报错并提示 `plink --make-bed`；混合 SNP+indel 模式（`6c 1b 00`，indel 占两个 2-bit hardcall）报错并提示去掉 indel。魔数共 **3 字节**，文件大小必须等于 `3 + 位点数 × ceil(样本数/4)`，否则报错——这条同时挡住 PLINK 在「单染色体位点数 > 524288」时改用的另一套布局。样本数非 4 的倍数时，尾部 padding 位忽略。
  - 2-bit 码：`00`=A1 纯合、`01`=缺失、`10`=杂合、`11`=A2 纯合。
- **`.bim`**：6 列 `chr / id / cm / bp / a1 / a2`；位置取 `bp` 列（`cm` 忽略）。等位若不是单个 `A/C/G/T`（indel、PLINK 的 `D`/`I`/`0` 记号）→ 该位点按 **Indel 过滤**丢弃（与 VCF 路径丢弃多字符 `REF`/`ALT` 同理）；`a1 == a2` → 按非双等位丢弃。
- **`.fam`**：6 列，**第 2 列 IID 为样本名**（与 `plink --recode vcf` 的 VCF 样本名一致），IID 重复 → 报错退出。
- **过滤与编码语义对齐 `-InVCF`（而非 `-InGenotype`）**：位点筛选规则、MAF/Miss/Het 分母、伪单倍型构造都与 VCF 路径一致，因此**位点集合与每个距离分箱的 `NumberPairs` 逐位相同**。但**在有相位的真实数据上 `r²`/`D′` 不与 `-InVCF` 逐字节一致**，原因是 2-bit `.bed` 不存相位：
  - 最佳/次佳等位由 `.bed` 计数决定（`00` 计作 GT 字符 `'0'`，`11` 计作 `'1'`，`10` 各计一次），并列时 `'0'` 为最佳——等价于 VCF 路径扫描 GT 字符。
  - `-Miss`/`-Het` 分母为**样本数**；`-MAF` 分母为**观测到的等位数**（`sed/(sed+max)`）。
  - 杂合（PLINK 无相位）→ 伪单倍型 `(best, sed)`，恒编码为 `(0,1)`；缺失 → `2`。
  - 与 `-InGenotype` 的已知差异（两条旧路径本就不同）：基因型路径按碱基字母选 best/sed、`-MAF` 分母为 `2×N`（含缺失）。故等位计数恰好相等的位点与 `-MAF` 边界位点，两条输入路径结果可能不同（golden 中记为 INFO）。
  - **相位（关键）**：`.bed` 只有基因型、没有相位，`-InPlink` 对每个位点按 `(best, sed)` 合成单倍型；`-InVCF` 用 `GT` 里的**真实相位**。1000G chr22 实测：位点集合与 300,000 个分箱的 `NumberPairs` 完全相同，但 `r²` 平均相差约 0.0055（最大 0.026）。等位朝向（PLINK 写 `A1=ALT`）不影响结果，因此这是相位效应而非等位顺序效应。**无相位数据（如合成的 `.geno`/随机配对输入）下两者逐字节一致**，这也是旧 golden 用例未暴露该差异的原因。
- **`-SubPop`**：样本名取 `.fam` IID（合成 `#CHROM POS <IID>...` 表头复用 `resolve_subpop_header`，`sample_start=2`）；多群体/两列写法与 VCF 一致。
- **`-L`**：与 VCF 回退路径相同——位点按所属区域分桶，跨区域对不计算。
- **`-OutFilterSNP`**：PLINK 无「原文行」概念，改为输出 `{prefix}.plink.filter.bim`（通过位点的 6 列文本，按 `.bim` 顺序，仅全样本路径）。因此该模式单线程读 `.bed`（保证顺序）。
- **`-EHH`**：与 VCF/基因型输入共用同一条 EHH 路径。
- **并行**：`.bed` 可随机访问，位点区间按 `-T` 切块并行解码（`≥4096` 位点才切，上限 64 块），各块缓冲按 worker 序拼接后仍与单线程的稳定排序/去重结果一致。
- **已知限制**：无索引快路径与「按位置分片」只对 `-InVCF` 生效（见 §10/§11），所以单染色体 PLINK 输入的 LD 计算与 `-InGenotype` 一样按 contig 串行（`-T` 仍并行化解码与多 contig 分配）。

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
- 旧 `Plot_OnePop.pl / Plot_MultiPop.pl` 无需修改即可读入。

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

## 输出：`.plink.filter.bim`（M7，`-InPlink -OutFilterSNP`）
- 通过全部过滤（Miss/Het/双等位/MAF）的 PLINK 位点，`chr / id / 0 / bp / a1 / a2` 6 列文本，**按 `.bim` 原顺序、无表头**；仅全样本路径产生（与其它输入的 SubPop 行为一致）。

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
- **每行以染色体名开头，行的长度不受缓冲区限制**：旧版用 C++ 流输出（无长度上限），本实现亦然（`src/rowfmt.h`）——放不下的行会改用堆缓冲重新格式化，而不是被截断。染色体名 ≳90 字符时旧实现会越界读并与旧版不一致，修复后逐字节一致。

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
- 数据行同样以染色体名开头，长度不受限制（见 `.LD.gz` 一节的说明）。
- 绘图失败会在 stderr 报一行 warning，但**退出码仍为 0**（与旧版一致）。