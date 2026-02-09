-------------------------
# biotools_summary
A repository that regularly collects and organizes bioinformatics-related tools.

## 目录
- [sylph](#sylph)：利用ANIs快速精准的物种水平宏基因组谱分析工具
- [GUS_Endobiotics](#gus_endobiotics)：识别肠道微生物β-葡糖醛酸酶
- [CRC_GUS](#crc_gus)：识别肠道微生物β-葡糖醛酸酶
- [Spacedust](#spacedust)：利用结构同源性与新型统计模型，实现微生物保守基因簇的De Novo发现
- [MMETHANE](#mmethane)：面向微生物组与代谢组的可解释 AI 预测模型
- [HMMER](#hmmer)：隐马尔可夫模型（HMM）序列分析工具套件，特别适用于蛋白质功能注释


--------------------------------------------------
## sylph
利用ANIs快速精准的物种水平宏基因组谱分析工具
宏基因组分析：鉴定样本中的物种/分类群及其丰度。
包含性平均核苷酸一致性查询：搜索一个基因组（例如大肠杆菌）是否存在于您的样本中，并估算其与查询基因组的ANI相似度。
其主要功能和特点包括：
1. 精确的物种级分析：比Kraken等工具假阳性更少，精度和灵敏度与MetaPhlAn等标记基因方法相当。
2. 极快的速度：可比其他方法快50倍以上，且支持多线程、多样本并行分析。
3. 内存效率高：分析整个GTDB-R220数据库（约11万个基因组）仅需约15GB内存。
4. 提供准确的ANI信息：即使基因组覆盖度低至0.1倍，也能给出准确的ANI估计值。
5. 数据库灵活：提供预构建的细菌、病毒、真核生物数据库，也支持用户轻松构建自定义数据库（如使用自己的MAGs）。
6. 读长兼容性好：同时适用于短读长和长读长测序数据，在牛津纳米孔的独立基准测试中表现最佳。

#### 文章引用
Shaw J, Yu YW. Rapid species-level metagenome profiling and containment estimation with sylph. Nat Biotechnol. 2025;43(8):1348-1359. doi:10.1038/s41587-024-02412-y

#### 文档手册
https://sylph-docs.github.io/

#### 实操命令
- 执行程序与路径
```sh
sylph=/home/data/t070605/Softwares/miniforge3/bin/sylph # 主程序
sylphtax=/home/data/t070605/Softwares/miniforge3/bin/sylph-tax # 物种id映射到metaphlan格式
db_sylph=/home/data/t070605/Databases/db_sylph # 预编译数据库
# or
sylph=/mnt/lustre/user/hanlj/Software/anaconda3/bin/sylph
sylphtax=/mnt/lustre/user/hanlj/Software/anaconda3/bin/sylph-tax
db_sylph=/mnt/lustre/user/hanlj/Databases/db_sylph
```

- sylph
```sh
# 可以指定一个或多个预编译数据库
sylph profile ${db_sylph}/gtdb-r226-c200-dbv1.syldb ${db_sylph}/imgvr_c200_v0.3.0.syldb -1 *_1.fastq.gz -2 *_2.fastq.gz -t 16 > profiling.tsv

# multi-sample single-end profiling
# sylph profile ${db_sylph}/gtdb-r226-c200-dbv1.syldb *.fastq -t 16 > profiling.tsv
```

- sylph-tax
```sh
# incorporate GTDB-r226 and IMGVR-4.1 taxonomies into sylph's results
sylph-tax taxprof sylph_results/*.tsv -t GTDB_r226 IMGVR_4.1 -o sylph_results/output_prefix-

# merge multiple results 
sylph-tax merge sylph_results/*.sylphmpa --column relative_abundance -o merged_abundance_file.tsv
```


--------------------------------------------------
## GUS_Endobiotics
识别肠道微生物β-葡糖醛酸酶，将输入的蛋白质FASTA序列与参考的细菌β-葡萄糖醛酸酶（GUS）序列进行比对，以评估每个输入的蛋白质序列是否满足现有残基位置和组成阈值，从而将其注释为推测的GUS蛋白。每个输入的蛋白质序列都会通过与其他已确认具有酶活性的细菌GUS酶参考序列进行比对，来评估在特定位置是否存在功能必需的氨基酸。

#### 文章引用
Simpson JB, Walker ME, Sekela JJ, et al. Gut microbial β-glucuronidases influence endobiotic homeostasis and are modulated by diverse therapeutics. Cell Host Microbe. 2024;32(6):925-944.e10. doi:10.1016/j.chom.2024.04.018

#### 文档手册
https://github.com/redinbolab/GUS_Endobiotics/tree/main
https://github.com/redinbolab/GUS_Endobiotics/tree/main/GUS%20Identification/GUS-ID

#### 实操命令（GUS Identification）
```sh
# 1. 安装必要的程序，依赖程序安装在checkm环境下
#dependencies:
#  - biopython=1.79.*
#  - blast=2.12.*
#  - python=3.11.*

# 2. 将GUS-ID_v2.py拷贝到分析目录
cp /home/data/t070605/Softwares/GUS_Endobiotics/GUS-ID_v2.py /path/to/analyse/


# 3. 依赖程序安装在checkm环境下，激活环境
conda activate checkm
/home/data/t070605/Softwares/miniforge3/envs/checkm/bin/python GUS-ID_v2.py {蛋白质输入文件.fasta}
```
脚本运行完毕后，将在当前目录生成以下结果文件：

- 汇总的 FASTA 文件：包含所有与参考序列匹配的 GUS 序列。
- 序列目录：为每个鉴定为 GUS 阳性的匹配序列生成独立的 FASTA 文件，以序列标题命名。
- 保守残基目录：包含详细说明查询序列与参考序列之间保守残基信息的 .txt 文件。
- 日志文件：记录脚本版本和运行时间等信息。

--------------------------------------------------
## CRC_GUS
肠道微生物β-葡萄糖醛酸酶在结直肠癌中作用的前沿研究。针对结直肠癌（CRC）中肠道微生物β-葡萄糖醛酸酶（gmGUSs）的作用，通过新型分析流程从临床样本中构建了包含550种酶的gmGUSs图谱。研究首次揭示了CRC发生发展中“微生物-gmGUS-代谢物”（MGM）轴的阶段特异性失调，该轴可通过重塑氨基酸代谢等途径促进肿瘤进展。分析发现，晚期CRC患者肠道中特定结构类型（如Mini-Loop 2）的gmGUSs以及来源于解纤维素拟杆菌（B. cellulosilyticus）的gmGUSs显著富集。基于差异gmGUSs构建的诊断模型（AUC 0.827）和预后评分模型（GUSscore）均表现出优异的预测性能，确立了gmGUSs作为潜在早期诊断标志物和预后评估指标的重要临床价值。

#### 文章引用
Chen J, Li Y, Tang S, Jin W, Yan R. Gut microbial β-glucuronidases and their role in the microbiome-metabolite axis in colorectal cancer. Nat Commun. 2025;16(1):10660. Published 2025 Nov 27. doi:10.1038/s41467-025-65679-y

#### 文档手册
https://github.com/yr2008-UM/CRC_GUS

#### 实操命令
1. GUS_Identification: Pipeline for GUS identification.
2. Loop_Classification: Pipeline for loop classification.
3. Taxonomic_Annotation: Pipeline for construction of kraken database.
4. DataAnalysis: Resources and R scripts for data analysis.

--------------------------------------------------
## Spacedust
利用结构同源性与新型统计模型，实现微生物保守基因簇的De Novo发现

#### Citing
Zhang R, Mirdita M, Söding J. De novo discovery of conserved gene clusters in microbial genomes with Spacedust. Nat Methods. 2025;22(10):2065-2073. doi:10.1038/s41592-025-02816-x

#### 文档手册
https://github.com/soedinglab/Spacedust/

#### 实操命令


--------------------------------------------------
## MMETHANE
面向微生物组与代谢组的可解释 AI 预测模型

#### 文章引用
Dawkins JJ, Gerber GK. MMETHANE: interpretable AI for predicting host status from microbial composition and metabolomics data. Microbiome. 2025;14(1):21. Published 2025 Dec 8. doi:10.1186/s40168-025-02270-z

#### 文档手册
https://github.com/gerberlab/mmethane

#### 实操命令



--------------------------------------------------
## HMMER
隐马尔可夫模型（HMM）序列分析工具套件，用于蛋白质序列同源性搜索和序列比对。HMMER 通过构建隐马尔可夫模型来捕捉蛋白质家族的多序列比对特征，能够更敏感地检测远缘同源序列。其中 `hmmscan` 命令用于将查询蛋白质序列与 HMM 数据库（如 Pfam）进行比对，特别适用于对宏基因组组装基因组（MAGs）或 contigs 的预测蛋白质进行功能注释。

用户经验：在 contigs 或 MAGs 序列经 prokka 预测基因后，直接使用 Pfam-A.hmm 数据库对预测的蛋白质序列进行同源基因检索，效果优于默认的 prokka 注释流程。

主要功能和特点包括：
1. **高灵敏度的同源检测**：利用隐马尔可夫模型捕捉蛋白质家族的保守模式，能够检测远缘同源序列。
2. **多种分析工具**：
   - `hmmbuild`：从多序列比对构建 HMM 模型
   - `hmmsearch`：使用 HMM 模型搜索蛋白质序列数据库
   - `hmmscan`：将查询序列与 HMM 数据库进行比对（适用于功能注释）
   - `hmmpress`：预处理 HMM 数据库以加速搜索
   - `jackhmmer`：迭代搜索，构建更全面的同源序列集
3. **高效的算法实现**：采用加速启发式算法（如 MSV、Bias 过滤器），在保持高灵敏度的同时大幅提升搜索速度。
4. **广泛的数据库兼容**：支持 Pfam、TIGRFAM、PANTHER 等常用 HMM 数据库。
5. **灵活的输入输出格式**：支持 FASTA、 Stockholm、 SELEX 等多种格式，输出结果可定制。

#### 文章引用
Eddy SR. Accelerated Profile HMM Searches. PLoS Comput Biol. 2011;7(10):e1002195. doi:10.1371/journal.pcbi.1002195

#### 文档手册
http://hmmer.org/
https://github.com/EddyRivasLab/hmmer

#### 实操命令
- **安装 HMMER**
```sh
# 通过 Conda 安装
conda install -c bioconda hmmer

# 或从源码编译
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.4
./configure
make
make install
```

- **下载 Pfam 数据库**
```sh
# 下载 Pfam-A.hmm 数据库
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# 预处理数据库（加速搜索）
hmmpress Pfam-A.hmm

# 数据库文件已下载并处理，路径如下：
ll /home/data/t070605/Databases/Pfam/
```

- **hmmscan 基本用法**
```sh
# 将预测的蛋白质序列与 Pfam 数据库比对
hmmscan --cpu 8 --domtblout output.domtblout Pfam-A.hmm protein_sequences.faa

# 常用参数说明：
# --cpu：使用的 CPU 线程数
# --domtblout：输出域表格结果（推荐格式）
# -E：E-value 阈值（默认 10.0）
# --domE：域 E-value 阈值
# --tblout：输出简洁的表格结果
# -o：标准输出文件
# --notextw：禁用文本换行，方便解析
```

- **结果解读**
`output.domtblout` 文件包含每个查询序列与 Pfam 家族的比对信息，主要列包括：
1. 查询序列名称
2. 目标 HMM 名称（Pfam ID）
3. E-value：统计显著性，值越小越显著
4. score：比对得分
5. bias：偏差校正值
6. 查询序列的匹配区间
7. HMM 模型的匹配区间

- **结合 prokka 输出的使用示例**
```sh
# prokka 预测基因后，提取蛋白质序列
# prokka 输出通常包含 .faa 文件（蛋白质 FASTA）

# 使用 hmmscan 进行功能注释
hmmscan --cut_ga --noali --cpu 1 --domtblout hmmscan_Blautia/hmmscan_L303A.txt /home/data/t070605/Databases/Pfam/Pfam-A.hmm ./03_prokka/L301A/L301A.faa
```
