# Spacedust：多基因组保守基因簇发现工具

Spacedust 是一个模块化工具包，用于基于同源性和基因邻接保守性，在多个基因组中识别保守基因簇。它整合了 Foldseek 的快速敏感结构比较能力和 MMseqs2 的同源搜索功能，引入了一种新颖的方法来聚合基因组对之间的同源命中集，并使用凝聚层次聚类算法识别具有保守基因邻接的命中簇。

**主要特性：**
- 支持基于序列（MMseqs2）和结构（Foldseek）的搜索模式
- 可处理核苷酸基因组（FASTA + GFF3）或蛋白质序列（FASTA）输入
- 支持迭代搜索提高灵敏度
- 可生成结构数据库映射到 AlphaFoldDB 或使用 ProstT5 预测 3D 结构
- 设计为在多核系统上高效运行
- 开源 GPLv3 许可，支持 Linux 和 macOS

## 安装示例

```bash
# 下载预编译版本（Linux AVX2）
wget https://mmseqs.com/spacedust/spacedust-linux-avx2.tar.gz
tar xvzf spacedust-linux-avx2.tar.gz
export PATH=$(pwd)/spacedust/bin/:$PATH

# 或使用 Conda 安装
conda install -c conda-forge -c bioconda spacedust

# 或从源码编译
git clone https://github.com/soedinglab/spacedust.git
cd spacedust
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make -j
make install
```

## 使用示例

### 1. 创建序列数据库

```bash
# 对于核苷酸输入（需 GFF3 注释）
spacedust createsetdb genome1.fna genome2.fna setDB tmpFolder \
  --gff-dir gffDir.txt --gff-type CDS

# 对于蛋白质输入
spacedust createsetdb genome1.faa genome2.faa setDB tmpFolder
```

### 2. 创建结构数据库

```bash
# 下载参考 Foldseek 数据库
foldseek databases Alphafold/UniProt refFoldseekDB tmpFolder

# 将序列数据库映射到结构数据库
spacedust aa2foldseek setDB refFoldseekDB tmpFolder
```

### 3. 运行聚类搜索

```bash
# 基本序列搜索（使用 MMseqs2）
spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder

# 结构+序列混合搜索
spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder \
  --search-mode 1

# 迭代搜索（类似 PSI-BLAST）
spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder \
  --num-iterations 2

# 仅使用 ProstT5 + Foldseek 搜索
spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder \
  --search-mode 2
```

### 4. 使用 ProstT5 创建结构数据库

```bash
# 下载 ProstT5 模型
foldseek databases ProstT5 weights tmpFolder

# 创建 Foldseek 数据库
foldseek createdb genome1.faa genome2.faa DB --prostt5-model weights

# 创建 Spacedust 数据库
spacedust createsetdb DB setDB tmpFolder
```

## 输出格式

Spacedust 输出制表符分隔的文本文件（.tsv），每个报告簇包含一个摘要行和多个命中行：

```
#clusterID  query_acc  target_acc  clusterMatchPvalue  multihitPvalue  num_hits
>queryID    targetID   bestHitPvalue  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen  alnCigar
```

摘要行以 `#` 开头，包含簇 ID、查询基因组、目标基因组、簇匹配 P 值、多重命中 P 值和命中数。每个后续行以 `>` 开头，描述簇中的单个成员命中。

## 依赖项

Spacedust 需要 Foldseek 进行结构比较。可将 Foldseek 二进制文件放在与 Spacedust 相同的目录中，或使用 `--foldseek-path` 参数指定路径。

## 参数说明

### `createsetdb` 命令参数
- `--gff-type`：GFF 文件中要过滤的特征类型（默认：""，所有特征）
- `--gff-dir`：GFF 目录文件路径

### `clustersearch` 命令参数
- `--search-mode`：0：使用 MMseqs2 进行序列搜索，1：使用 Foldseek 进行结构比较，2：Foldseek + ProstT5（默认：0）
- `--num-iterations`：迭代轮廓搜索迭代次数（默认：1）
- `--profile-cluster-search`：执行轮廓（目标）-序列搜索
- `--filter-self-match`：移除同一集合间的命中
- `--max-gene-gap`：两个簇合并之间允许的最大基因数（默认：3）
- `--cluster-size`：定义簇的最小基因数（默认：2）
- `--foldseek-path`：Foldseek 二进制文件路径

## 临时文件管理

在执行过程中，Spacedust 会将所有中间输出保存在 `tmpFolder` 中。可以通过传递 `--remove-tmp-files` 参数在流程完成后清理临时文件夹。

---

**参考文献**：Zhang, R., Mirdita, M., & Söding, J. (2024). De novo discovery of conserved gene clusters in microbial genomes with Spacedust. bioRxiv. doi:10.1101/2024.10.02.616292

**GitHub 仓库**：https://github.com/soedinglab/Spacedust

**许可证**：GPLv3

**语言**：C++

**最后更新**：2026年1月30日