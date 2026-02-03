-------------------------
# biotools_summary
A repository that regularly collects and organizes bioinformatics-related tools.

## 目录
- [sylph](#sylph)：利用ANIs快速精准的物种水平宏基因组谱分析工具
- [GUS_Endobiotics](#gus_endobiotics)：识别肠道微生物β-葡糖醛酸酶
- [Spacedust](#spacedust)：利用结构同源性与新型统计模型，实现微生物保守基因簇的De Novo发现
- [MMETHANE](#mmethane)：面向微生物组与代谢组的可解释 AI 预测模型


--------------------------------------------------
## sylph
利用ANIs快速精准的物种水平宏基因组谱分析工具

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
识别肠道微生物β-葡糖醛酸酶

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

#### 文章引用

#### 文档手册
https://github.com/yr2008-UM/CRC_GUS

#### 实操命令


--------------------------------------------------
## Spacedust
利用结构同源性与新型统计模型，实现微生物保守基因簇的De Novo发现

#### Citing
Zhang R, Mirdita M, Söding J. De novo discovery of conserved gene clusters in microbial genomes with Spacedust. Nat Methods. 2025;22(10):2065-2073. doi:10.1038/s41592-025-02816-x

#### 文档手册
https://github.com/soedinglab/Spacedust/

#### 实操命令
01.GUS_Identification: Pipeline for GUS identification.
02.Loop_Classification: Pipeline for loop classification.
03.Taxonomic_Annotation: Pipeline for construction of kraken database.
04.DataAnalysis: Resources and R scripts for data analysis.

--------------------------------------------------
## MMETHANE
面向微生物组与代谢组的可解释 AI 预测模型

#### 文章引用
Dawkins JJ, Gerber GK. MMETHANE: interpretable AI for predicting host status from microbial composition and metabolomics data. Microbiome. 2025;14(1):21. Published 2025 Dec 8. doi:10.1186/s40168-025-02270-z

#### 文档手册
https://github.com/gerberlab/mmethane

#### 实操命令



