# GUS_Endobiotics：肠道微生物β-葡萄糖醛酸苷酶

本仓库包含与**《肠道微生物β-葡萄糖醛酸苷酶影响内源性物质稳态并受多种治疗药物调节》**相关的所有代码，该文章已被《细胞·宿主与微生物》接受发表。

## 概述

GUS_Endobiotics 是一个用于识别和分析肠道微生物β-葡萄糖醛酸苷酶的综合工具包。它提供以下工具：
- 从蛋白质FASTA文件中识别GUS酶
- 按结构环类型对GUS蛋白进行分类
- 分析结构信息指导的速率-丰度关联
- 执行速率-微生物丰度关联的统计分析

## 仓库结构

```
GUS_Endobiotics/
├── GUS Identification/                # GUS识别
│   ├── GUS-ID/                        # GUS识别流程
│   │   ├── GUS-ID_v1.py               # GUS识别脚本版本1
│   │   ├── GUS-ID_v2.py               # GUS识别脚本版本2
│   │   ├── GUS_test.fa                # 测试蛋白质序列
│   │   ├── create_GUS-ID_v2_env.yml   # GUS-ID_v2的Conda环境
│   │   ├── open_to_break_pkg/         # 包含参考序列的包
│   │   └── GUS-ID_readme_v2.txt       # 详细使用说明
│   └── GUS Structural Class ID/       # GUS环分类
│       ├── GUS_Loop_ID_v3.2.py        # 环分类脚本
│       └── GUS_Reference_Sequences/   # 分类参考序列
├── Proteomic Abundance-Rate Regression Analysis/  # FitFinder工具
├── Omics-Taxonomic Rank t-Tests/                  # 统计分析
└── Various figure generation scripts/             # 论文图表生成脚本
```

## 安装

### 先决条件
- **Miniconda** >= 4.10.3 (推荐)
- **Python** >= 3.11
- **BLAST** >= 2.12.0
- **BioPython** >= 1.79.0

### 设置GUS-ID环境

1. **标准架构：**
```bash
conda env create -n GUS-ID_v2 --file create_GUS-ID_v2_env.yml
conda activate GUS-ID_v2
```

2. **macOS arm64架构：**
```bash
# 单独下载适用于osx-arm64的BLASTp 2.12.0
curl -L https://anaconda.org/bioconda/blast/2.12.0/download/osx-64/blast-2.12.0-h0370960_3.tar.bz2 -o blast-2.12.0.tar.bz2

# 创建环境
conda env create -n GUS-ID_v2 --file create_GUS-ID_v2_env_arm64.yml
conda activate GUS-ID_v2

# 安装BLASTp
conda install blast-2.12.0.tar.bz2
```

### 环境依赖 (create_GUS-ID_v2_env.yml)
```yaml
name: GUS-ID_v2.0
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - biopython=1.79.*
  - blast=2.12.*
  - python=3.11.*
```

## 使用方法

### 1. GUS识别 (GUS-ID)

**目的：** 通过将蛋白质FASTA序列与参考细菌GUS序列进行比较，并评估功能必需氨基酸的保守性，从而识别β-葡萄糖醛酸苷酶。

**输入要求：**
- 蛋白质FASTA文件 (`.fa` 或 `.fasta` 扩展名)
- 每个序列必须具有唯一的标题
- 无大小限制，但对于大于几百MB的文件，请使用GNU parallel

**基本用法：**
```bash
# 激活环境
conda activate GUS-ID_v2

# 运行GUS识别
python GUS-ID_v2.py input_proteins.fasta

# 使用提供的示例测试
python GUS-ID_v2.py GUS_test.fa
```

**输出：**
- `{输入名称}_GUS_aligned/` - 比对文件
- `{输入名称}_GUS_conserved_residues/` - 详细描述保守残基的文本文件
- `{输入名称}_GUS_individual_seqs/` - 每个GUS阳性匹配的单独FASTA文件
- 所有GUS匹配序列的拼接FASTA文件
- 包含版本和运行时间信息的日志文件

### 2. GUS结构环分类 (GUS_Loop_ID)

**目的：** 基于与参考GUS序列的多序列比对，对GUS蛋白按其结构环类型进行分类。

**用法：**
```bash
python GUS_Loop_ID_v3.2.py input_proteins.fasta
```

**脚本概述：**
```python
# GUS_Loop_ID_v3.2.py中的关键参数
in_file = sys.argv[1]  # 输入FASTA文件
ofnm_afa = in_file.replace('.fasta','').replace('.fa','')+'.afa'  # 输出比对文件
csv_filename = in_file.replace('.fasta','').replace('.fa','') + '_Loop_Classifications.csv'

# 主要工作流程：
# 1. 将输入序列与参考序列拼接
# 2. 使用ClustalOmega进行多序列比对
# 3. 基于比对位置对环进行分类
# 4. 输出包含环分类的CSV文件
```

### 3. GUS-ID脚本代码结构

**GUS-ID_v2.py的主要组件：**
```python
#!/usr/bin/env python
import sys
import os
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Blast import NCBIXML
import open_to_break_pkg.util

# 比对阈值 (25% 同一性)
min_alignment_identity_threshold = 0.25

# 输入/输出配置
seed_fasta_dir = 'open_to_break_pkg/GUS_reference_sequences'
db_fasta_file_str = sys.argv[1]  # 输入FASTA文件

# 输出目录
alignment_root_out_path = db_fasta_file_name + '_GUS_aligned'
residues_root_out_path = db_fasta_file_name + '_GUS_conserved_residues'
fasta_root_out_path = db_fasta_file_name + '_GUS_individual_seqs'

# 核心功能：将查询序列与参考GUS蛋白比对
def align_all_seeds(seed_id2proteins, candidate, min_alignment_identity_threshold,
                    alignment_root_out_path, residues_root_out_path, fasta_root_out_path):
    # 执行BLAST比对和残基保守性检查
    # 如果所有7个保守位置都存在，则返回匹配残基
```

## 示例工作流程

### 完整的GUS分析流程

```bash
# 步骤1：设置环境
conda create -n gus_analysis python=3.11 biopython=1.79 blast=2.12
conda activate gus_analysis

# 步骤2：从宏基因组蛋白质中识别GUS酶
python GUS-ID_v2.py metagenome_proteins.fasta

# 步骤3：按环类型对识别的GUS蛋白进行分类
python GUS_Loop_ID_v3.2.py metagenome_proteins_GUS_individual_seqs/*.fasta

# 步骤4：分析结果
# - 在 *_GUS_conserved_residues/ 目录中检查保守残基
# - 在 *_Loop_Classifications.csv 中查看环分类
# - 使用单个序列进行下游分析
```

## 主要特性

1. **全面的GUS识别**：使用BLAST比对，阈值为25%同一性，并检查7个关键残基位置的保守性。

2. **结构分类**：根据结构特征将GUS蛋白分类为环类型（无环、迷你环、环1、环2等）。

3. **质量控制**：
   - 验证输入序列长度相对于参考序列的长度
   - 检查序列标题的唯一性
   - 提供详细的残基保守性报告

4. **灵活的输入**：适用于任何蛋白质FASTA文件，从单个序列到大型宏基因组目录。

5. **可重现的环境**：Conda环境文件确保跨平台依赖关系的一致性。

## 输出解读

### GUS识别输出：
- **保守残基文件**：文本文件，显示每个查询-参考对的比对位置和残基匹配情况。
- **单个序列文件**：每个确认GUS的FASTA文件，可用于系统发育分析。
- **比对文件**：BLAST XML格式的比对文件（处理过程中临时生成）。

### 环分类输出：
- **CSV文件**：列包括蛋白质ID、环类别（Loop1, Loop2）和总环类别。
- **比对文件**：ClustalOmega多序列比对，AFA格式。

## 高级用法

### 大型数据集的并行处理
对于超过几百MB的FASTA文件，使用GNU parallel：
```bash
# 拆分输入文件并并行处理
python split_fasta.py large_input.fasta 100  # 拆分为100个序列的块
parallel -j 8 python GUS-ID_v2.py {} ::: chunk_*.fasta
```

### 自定义参考数据库
要使用自定义GUS参考序列：
1. 替换 `open_to_break_pkg/GUS_reference_sequences/` 中的序列
2. 在 `open_to_break_pkg/GUS_seed_indices.txt` 中更新新的残基位置
3. 使用更新后的参考序列重新运行GUS-ID

## 故障排除

### 常见问题：

1. **"BLAST not found" 错误**：确保BLAST已安装在conda环境中：
   ```bash
   conda activate GUS-ID_v2
   which blastp  # 应显示conda环境中的路径
   ```

2. **大文件内存问题**：使用序列过滤或拆分输入文件。

3. **重复的序列标题**：确保所有输入序列具有唯一的ID。

4. **macOS arm64兼容性**：遵循BLAST安装的特殊说明。

## 引用

如果在出版物中使用此软件，请引用：

Simpson JB, Walker ME, Sekela JJ, 等. Gut microbial β-glucuronidases influence endobiotic homeostasis and are modulated by diverse therapeutics. *Cell Host & Microbe*. 2024;32(6):925-944.e10. doi:10.1016/j.chom.2024.04.018

## 许可证与联系

**重要提示**：此代码完全归北卡罗来纳大学教堂山分校Matthew Redinbo实验室所有 (https://www.redinbolab.org/)。

**预期用途**：此代码仅用于学术研究。

**联系**：对于运行此代码的问题或关于出版物使用的疑问，请联系 redinbo@unc.edu。

**引用咨询**：如果打算发表使用此代码收集的数据，请咨询 redinbo@unc.edu 以获取正确的引用方式。

## 仓库信息

- **GitHub**: https://github.com/redinbolab/GUS_Endobiotics
- **主分支**: https://github.com/redinbolab/GUS_Endobiotics/tree/main
- **GUS-ID文档**: https://github.com/redinbolab/GUS_Endobiotics/tree/main/GUS%20Identification/GUS-ID
- **最后更新**: 基于截至2025年的仓库状态
