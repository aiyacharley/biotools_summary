# HMMER：隐马尔可夫模型序列分析工具套件

HMMER 是一个基于隐马尔可夫模型（Hidden Markov Model, HMM）的生物信息学工具套件，主要用于蛋白质序列的同源性搜索、多序列比对和蛋白质家族分析。该工具套件由 Sean Eddy 实验室开发，已成为蛋白质功能注释、结构域识别和系统发育分析的标准工具之一。

**用户经验分享**：在对 contigs 或宏基因组组装基因组（MAGs）进行基因预测后（如使用 prokka），直接使用 `hmmscan` 命令配合 Pfam-A.hmm 数据库对预测的蛋白质序列进行功能注释，效果显著优于默认的 prokka 注释流程，能够识别更多远缘同源基因。

**主要应用场景**：
- 蛋白质序列的功能注释
- 蛋白质结构域识别
- 蛋白质家族成员鉴定
- 远缘同源序列检测
- 宏基因组/宏转录组数据分析
- 比较基因组学研究

## 主要特性

### 1. 高灵敏度的同源检测
- 利用隐马尔可夫模型捕捉蛋白质家族的保守模式
- 能够检测远缘同源序列（序列相似性低至 20-30%）
- 比传统 BLAST 方法更敏感，特别是对于结构域水平的相似性

### 2. 完整的工具套件
- **`hmmbuild`**：从多序列比对构建 HMM 模型
- **`hmmsearch`**：使用 HMM 模型搜索蛋白质序列数据库
- **`hmmscan`**：将查询序列与 HMM 数据库进行比对（功能注释）
- **`hmmpress`**：预处理 HMM 数据库以加速搜索
- **`jackhmmer`**：迭代搜索，构建更全面的同源序列集
- **`hmmalign`**：将序列与 HMM 模型进行比对
- **`hmmconvert`**：HMM 格式转换工具

### 3. 高效的算法实现
- 采用加速启发式算法（MSV、Bias 过滤器）
- 在保持高灵敏度的同时大幅提升搜索速度
- 支持多线程并行计算
- 内存效率优化，适合大规模数据分析

### 4. 广泛的数据库兼容
- **Pfam**：蛋白质家族数据库（最常用）
- **TIGRFAM**：蛋白质功能家族数据库
- **PANTHER**：蛋白质分类和功能预测数据库
- **SUPERFAMILY**：结构域超家族数据库
- **自定义数据库**：用户可构建特定项目的 HMM 数据库

### 5. 灵活的输入输出格式
- **输入格式**：FASTA、Stockholm、SELEX、Clustal、A2M
- **输出格式**：表格输出、域表格输出、多序列比对格式
- **结果可定制**：支持多种统计阈值和过滤选项

## 系统要求

### 基本要求
- **操作系统**：Linux、macOS、Windows（通过 WSL 或 Cygwin）
- **内存**：建议 4GB+（大型数据库如 Pfam 需要更多内存）
- **存储**：Pfam 数据库约 2-3GB（压缩后）
- **CPU**：支持多线程，建议多核处理器

### 依赖软件
- **C 编译器**（如 gcc）：用于源码编译
- **zlib 开发库**：数据压缩支持
- **OpenMP**：多线程支持（可选但推荐）
- **Perl/Python**：用于结果处理和脚本编写

## 安装指南

### 方法一：Conda 安装（推荐）
```bash
# 创建并激活环境
conda create -n hmmer -c bioconda hmmer
conda activate hmmer

# 验证安装
hmmscan -h
```

### 方法二：源码编译安装
```bash
# 下载最新版本（当前为 HMMER 3.4）
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar zxf hmmer-3.4.tar.gz
cd hmmer-3.4

# 配置和编译
./configure --prefix=/your/install/path
make
make install

# 添加环境变量
export PATH=/your/install/path/bin:$PATH
```

### 方法三：包管理器安装
```bash
# Ubuntu/Debian
sudo apt-get install hmmer

# CentOS/RHEL
sudo yum install hmmer

# macOS (Homebrew)
brew install hmmer
```

## 数据库准备

### Pfam 数据库下载和预处理
```bash
# 下载 Pfam-A.hmm 数据库（最新版本）
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz

# 解压缩
gunzip Pfam-A.hmm.gz

# 预处理数据库（加速搜索，必需步骤）
hmmpress Pfam-A.hmm

# 预处理生成的文件：
# Pfam-A.hmm.h3m - 主 HMM 文件
# Pfam-A.hmm.h3i - HMM 索引
# Pfam-A.hmm.h3f - HMM 数据文件
# Pfam-A.hmm.h3p - HMM 配置文件
```

### 其他常用数据库
```bash
# TIGRFAM 数据库
wget https://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz

# PANTHER 数据库
wget http://data.pantherdb.org/ftp/hmm_library/current_release/PANTHER.hmm.gz

# 解压和预处理
gunzip TIGRFAMs_15.0_HMM.LIB.gz
hmmpress TIGRFAMs_15.0_HMM.LIB
```

## 主要命令详解

### 1. `hmmscan` - 序列对 HMM 数据库搜索
**主要用途**：将查询蛋白质序列与 HMM 数据库比对，用于功能注释。

```bash
# 基本用法
hmmscan [选项] <hmm数据库> <查询序列文件>

# 常用参数
--cpu <n>                # 使用的 CPU 线程数
-E <float>              # E-value 阈值（默认 10.0）
--domE <float>          # 域 E-value 阈值
--incE <float>          # 包含性 E-value 阈值
--incdomE <float>       # 包含性域 E-value 阈值
-T <float>              # 比特分数阈值
--domT <float>          # 域比特分数阈值
--tblout <file>         # 输出简洁表格
--domtblout <file>      # 输出域表格（推荐）
--pfamtblout <file>     # 输出 Pfam 格式表格
-o <file>               # 标准输出文件
--notextw               # 禁用文本换行，方便解析
--cut_ga                # 使用 GA（收集）阈值
--cut_nc                # 使用 NC（噪声）阈值
--cut_tc                # 使用 TC（可信）阈值
```

### 2. `hmmsearch` - HMM 对序列数据库搜索
**主要用途**：使用已知的 HMM 模型搜索蛋白质序列数据库。

```bash
# 基本用法
hmmsearch [选项] <hmm模型> <序列数据库>

# 常用参数（与 hmmscan 类似）
--cpu <n>                # CPU 线程数
-E <float>              # E-value 阈值
--tblout <file>         # 表格输出
--domtblout <file>      # 域表格输出
```

### 3. `hmmbuild` - 构建 HMM 模型
**主要用途**：从多序列比对构建新的 HMM 模型。

```bash
# 基本用法
hmmbuild [选项] <输出hmm文件> <输入比对文件>

# 常用参数
--amino                # 输入为氨基酸序列（默认）
--dna                  # 输入为 DNA 序列
--rna                  # 输入为 RNA 序列
--fast                 # 快速模式（降低精度）
--symfrac <float>      # 符号分数阈值（默认 0.5）
--fragthresh <float>   # 片段阈值（默认 0.5）
--wpb                  # 加权位置背景
--wgsc                 # 加权全局序列权重
--wid                  # 加权 identity 权重
--eent                 # 有效熵权重
--eset                 # 有效集权重
--ere                  # 有效相对熵权重
--popen <float>        # 开放概率（默认 0.02）
--pextend <float>      # 扩展概率（默认 0.4）
```

### 4. `hmmpress` - 预处理 HMM 数据库
**主要用途**：创建 HMM 数据库的二进制索引，加速搜索。

```bash
# 基本用法
hmmpress <hmm数据库文件>

# 生成文件
<hmmdb>.h3m           # HMM 主文件
<hmmdb>.h3i           # HMM 索引文件
<hmmdb>.h3f           # HMM 数据文件
<hmmdb>.h3p           # HMM 配置文件
```

### 5. `jackhmmer` - 迭代搜索
**主要用途**：通过迭代搜索构建更全面的同源序列集。

```bash
# 基本用法
jackhmmer [选项] <查询序列> <目标数据库>

# 常用参数
-N <n>                 # 迭代次数（默认 1）
-E <float>             # 每轮迭代的 E-value 阈值
--incE <float>         # 包含性 E-value 阈值
--F1 <float>           # 第一轮序列分数阈值
--F2 <float>           # 第二轮序列分数阈值
--F3 <float>           # 第三轮序列分数阈值
```

## 实战案例

### 案例一：MAGs/Contigs 蛋白质功能注释
```bash
# 1. 从 prokka 输出提取蛋白质序列
# prokka 输出目录包含 .faa 文件（蛋白质 FASTA）

# 2. 使用 hmmscan 进行功能注释
hmmscan \
  --cpu 16 \
  --domtblout mags_proteins.domtblout \
  --tblout mags_proteins.tblout \
  -E 1e-5 \
  --notextw \
  Pfam-A.hmm \
  prokka_output/*.faa

# 3. 提取显著匹配（E-value < 1e-5）
awk '!/^#/ && $7 < 1e-5' mags_proteins.domtblout > significant_hits.tsv

# 4. 按 Pfam 家族统计
cut -f2 significant_hits.tsv | sort | uniq -c | sort -nr > pfam_counts.txt

# 5. 生成每个查询序列的最佳匹配
awk '!/^#/ {print $1 "\t" $2 "\t" $7 "\t" $4}' mags_proteins.domtblout | \
  sort -k1,1 -k3,3g | \
  awk '!seen[$1]++' > best_hits_per_gene.tsv
```

### 案例二：特定蛋白质家族的鉴定
```bash
# 1. 构建特定家族的 HMM 模型（如需）
# 从 Pfam 下载特定家族或自己构建

# 2. 搜索特定家族成员
hmmsearch \
  --cpu 8 \
  --domtblout ABC_transporters.domtblout \
  -E 1e-10 \
  ABC_transporter.hmm \
  all_proteins.faa

# 3. 提取匹配序列
grep -v '^#' ABC_transporters.domtblout | cut -f1 | sort -u > abc_transporter_ids.txt
seqtk subseq all_proteins.faa abc_transporter_ids.txt > abc_transporter_seqs.faa
```

### 案例三：宏转录组数据分析
```bash
# 1. 转录组组装和翻译
# Trinity 组装后，TransDecoder 预测蛋白质

# 2. 功能注释
hmmscan \
  --cpu 32 \
  --domtblout metatranscriptome.domtblout \
  -E 1e-3 \
  --incE 0.01 \
  Pfam-A.hmm \
  transdecoder_output/longest_orfs.pep

# 3. 功能分类统计
# 使用 Pfam 分类信息进行功能归类
```

### 案例四：自定义数据库构建
```bash
# 1. 收集同源序列
# 从 UniProt、NCBI 等下载相关序列

# 2. 多序列比对
mafft --auto input_sequences.fasta > aligned.fasta

# 3. 构建 HMM 模型
hmmbuild --amino custom_family.hmm aligned.fasta

# 4. 校准模型（可选）
hmmcalibrate custom_family.hmm

# 5. 使用自定义模型搜索
hmmsearch --cpu 8 --tblout results.tblout custom_family.hmm target_proteins.faa
```

## 输入输出格式详解

### 输入格式
1. **FASTA 格式**：标准的蛋白质序列格式
   ```
   >sequence_id description
   MKLLKTFLLSLFSLFSLCCLLLSCLCCLLLSCLCCLLLSCLCCLLLSCL
   ```
2. **Stockholm 格式**：多序列比对格式，包含序列和注释
3. **SELEX 格式**：另一种多序列比对格式
4. **HMM 格式**：HMMER 专用格式，包含模型参数

### 输出格式
1. **标准输出**（-o）：详细的文本输出，包含比对信息
2. **表格输出**（--tblout）：简洁的表格格式，包含基本统计
3. **域表格输出**（--domtblout）：**推荐格式**，包含域级别信息

#### domtblout 文件格式解析
```
# 列说明：
1. target_name         # 目标 HMM 名称（如 PF00001）
2. target_accession    # 目标 HMM 登录号
3. query_name          # 查询序列名称
4. query_accession     # 查询序列登录号
5. full_sequence_E-value  # 全序列 E-value
6. full_sequence_score    # 全序列比特分数
7. full_sequence_bias     # 全序列偏差
8. domain_number       # 域编号（同一序列可能有多个域）
9. domain_total        # 该序列匹配的总域数
10. domain_c-Evalue    # 条件 E-value（域水平）
11. domain_i-Evalue    # 独立 E-value（域水平）
12. domain_score       # 域比特分数
13. domain_bias        # 域偏差
14. hmm_from           # HMM 匹配起始位置
15. hmm_to             # HMM 匹配结束位置
16. ali_from           # 比对起始位置
17. ali_to             # 比对结束位置
18. env_from           # 包络起始位置
19. env_to             # 包络结束位置
20. acc                # 比对精度
21. target_description # 目标描述
```

## 结果解读和过滤

### 统计显著性指标
1. **E-value**：期望值，表示随机匹配的期望数量。值越小越显著。
   - `< 1e-10`：非常显著
   - `< 1e-5`：显著
   - `< 0.01`：可能显著
   - `> 0.01`：需要谨慎对待

2. **比特分数（bit score）**：比对质量的度量，与序列长度无关。
   - 越高表示比对质量越好
   - 通常与 E-value 结合使用

3. **条件 E-value（c-Evalue）**：考虑同一次搜索中其他序列的 E-value。
4. **独立 E-value（i-Evalue）**：不考虑其他序列的 E-value。

### 过滤策略
```bash
# 1. 基于 E-value 过滤
awk '!/^#/ && $7 < 1e-5' results.domtblout > filtered_results.tsv

# 2. 基于比特分数过滤
awk '!/^#/ && $6 > 50' results.domtblout > high_score_hits.tsv

# 3. 提取每个查询序列的最佳匹配
awk '!/^#/ {print $1 "\t" $3 "\t" $7 "\t" $6}' results.domtblout | \
  sort -k2,2 -k3,3g | \
  awk '!seen[$2]++' > best_hits.tsv

# 4. 去除重复匹配（同一序列匹配同一家族）
awk '!seen[$1,$3]++' results.domtblout > unique_hits.tsv
```

## 性能优化

### 1. 数据库预处理
```bash
# 预处理可显著加速搜索
hmmpress Pfam-A.hmm

# 预处理后搜索速度提升 2-5 倍
```

### 2. 并行处理
```bash
# 使用 GNU parallel 处理多个样本
parallel -j 4 "hmmscan --cpu 4 --domtblout {}.domtblout Pfam-A.hmm {}.faa" ::: sample1 sample2 sample3 sample4

# 或使用批处理脚本
for sample in *.faa; do
  base=$(basename $sample .faa)
  hmmscan --cpu 8 --domtblout ${base}.domtblout Pfam-A.hmm $sample &
done
wait
```

### 3. 参数调优
```bash
# 调整 E-value 阈值
hmmscan -E 1e-3 --incE 0.01 ...  # 宽松过滤
hmmscan -E 1e-10 --incE 1e-5 ... # 严格过滤

# 使用 GA/TC/NC 阈值（数据库提供时）
hmmscan --cut_ga ...  # 使用收集阈值
hmmscan --cut_tc ...  # 使用可信阈值
hmmscan --cut_nc ...  # 使用噪声阈值
```

### 4. 内存和磁盘优化
```bash
# 使用临时文件处理大数据库
hmmscan --tmpfile /tmp/hmmer_tmp ...

# 限制内存使用
hmmscan --max ...

# 分批处理大文件
split -l 1000 large_proteins.faa proteins_chunk_
for chunk in proteins_chunk_*; do
  hmmscan --cpu 4 --domtblout ${chunk}.domtblout Pfam-A.hmm $chunk
done
```

## 与其他工具的整合

### 1. 与 prokka 整合
```bash
# 完整的工作流程
# 1. 基因预测
prokka --outdir prokka_output --prefix mag contigs.fasta

# 2. 提取蛋白质序列
cp prokka_output/mag.faa .

# 3. HMMER 功能注释
hmmscan --cpu 16 --domtblout mag.domtblout Pfam-A.hmm mag.faa

# 4. 生成注释表格
awk '!/^#/ {print $3 "\t" $1 "\t" $21}' mag.domtblout | sort -u > mag_annotations.tsv
```

### 2. 与 eggNOG-mapper 结合
```bash
# 1. eggNOG-mapper 初步注释
emapper.py -i proteins.faa -o eggnog_output --cpu 8

# 2. HMMER 补充注释（特定家族）
hmmscan --cpu 8 --domtblout specific_families.domtblout custom_families.hmm proteins.faa

# 3. 结果整合
```

### 3. 与 InterProScan 对比
```bash
# InterProScan 是综合注释工具
interproscan.sh -i proteins.faa -f tsv -o ipr_output.tsv

# HMMER 可用于特定分析或自定义数据库
```

## 常见问题解决

### 1. 内存不足
```
错误：out of memory
解决方案：
- 使用分批处理
- 增加系统内存
- 使用 --max 参数限制内存
```

### 2. 搜索速度慢
```
问题：大型数据库搜索耗时
解决方案：
- 确保数据库已预处理（hmmpress）
- 增加 CPU 线程数（--cpu）
- 使用更严格的 E-value 阈值
- 考虑使用硬件加速（如有）
```

### 3. 结果太多/太少
```
调整策略：
- 太多结果：降低 E-value 阈值（如 1e-10）
- 太少结果：提高 E-value 阈值（如 1e-3）
- 使用 --incE 控制包含性阈值
```

### 4. 数据库格式问题
```
错误：invalid HMM file format
解决方案：
- 确保数据库文件完整下载
- 检查文件是否为 gzip 压缩（需要解压）
- 重新下载数据库文件
```

### 5. 多序列比对构建 HMM 失败
```
问题：hmmbuild 产生警告
解决方案：
- 检查输入比对文件格式
- 确保序列已正确比对
- 尝试不同的权重方案（--wpb, --wgsc 等）
```

## 进阶应用

### 1. 宏基因组 bin 质量评估
```bash
# 使用保守标记基因评估 MAG 完整性
# 下载单拷贝基因 HMM 数据库（如 BUSCO、CheckM）
hmmscan --cpu 8 --domtblout scg.domtblout single_copy_genes.hmm mag_proteins.faa

# 统计单拷贝基因存在情况
grep -c ">" single_copy_genes.hmm  # 总基因数
grep -v "^#" scg.domtblout | cut -f1 | sort -u | wc -l  # 检测到的基因数
```

### 2. 水平基因转移检测
```bash
# 通过非同源基因分布检测 HGT
# 1. 鉴定核心基因和附属基因
# 2. 分析基因在不同物种中的分布
# 3. 使用 HMMER 进行同源基因鉴定
```

### 3. 代谢通路重构
```bash
# 使用 KEGG Orthology (KO) 数据库
# 1. 下载 KO HMM 数据库
# 2. 鉴定代谢相关基因
# 3. 重构代谢通路
```

## 参考文献

1. **核心文献**：
   - Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol*. 2011;7(10):e1002195. doi:10.1371/journal.pcbi.1002195
   - Eddy SR. Profile hidden Markov models. *Bioinformatics*. 1998;14(9):755-763. doi:10.1093/bioinformatics/14.9.755

2. **应用文献**：
   - Finn RD et al. The Pfam protein families database: towards a more sustainable future. *Nucleic Acids Res*. 2016;44(D1):D279-D285. doi:10.1093/nar/gkv1344
   - Haft DH et al. TIGRFAMs: a protein family resource for the functional identification of proteins. *Nucleic Acids Res*. 2003;31(1):371-373. doi:10.1093/nar/gkg128

3. **数据库资源**：
   - **Pfam**: http://pfam.xfam.org/
   - **TIGRFAM**: https://www.jcvi.org/research/tigrfams
   - **PANTHER**: http://www.pantherdb.org/
   - **SUPERFAMILY**: http://supfam.org/

## 许可证和引用

### 许可证
HMMER 在 GNU General Public License (GPL) 版本 3 下发布。

### 引用 HMMER
使用 HMMER 时请引用：
```
Eddy SR. Accelerated Profile HMM Searches. PLoS Comput Biol. 2011;7(10):e1002195.
```

### 官方资源
- **官方网站**: http://hmmer.org/
- **GitHub 仓库**: https://github.com/EddyRivasLab/hmmer
- **文档**: http://eddylab.org/software/hmmer/Userguide.pdf
- **邮件列表**: hmmer-announce@ebi.ac.uk

## 更新日志

### 版本信息
- **当前稳定版**: HMMER 3.4（2021年发布）
- **主要改进**: 性能优化、新算法、更好的多线程支持

### 项目维护
- **维护者**: Sean Eddy 实验室
- **活跃开发**: 是
- **社区支持**: 活跃的邮件列表和 GitHub 问题跟踪

---

**最后更新**: 2026年2月9日
**文档作者**: 根据用户经验整理
**适用版本**: HMMER 3.x
**关键词**: 隐马尔可夫模型、蛋白质功能注释、hmmscan、Pfam、宏基因组分析