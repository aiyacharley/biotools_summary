# Sylph - fast and precise species-level metagenomic profiling with ANIs

> 所有关于 sylph 的文档已经迁移至 https://sylph-docs.github.io/.

## Sylph 是一款用于超快速、精确宏基因组分析的软件。
宏基因组分析：鉴定样本中的物种/分类群及其丰度。
包含性平均核苷酸一致性查询：搜索一个基因组（例如大肠杆菌）是否存在于您的样本中，并估算其与查询基因组的ANI相似度。

#### 其主要功能和特点包括：

1. 精确的物种级分析：比Kraken等工具假阳性更少，精度和灵敏度与MetaPhlAn等标记基因方法相当。

2. 极快的速度：可比其他方法快50倍以上，且支持多线程、多样本并行分析。

3. 内存效率高：分析整个GTDB-R220数据库（约11万个基因组）仅需约15GB内存。现在已更新至GTDB-R226版本。

4. 提供准确的ANI信息：即使基因组覆盖度低至0.1倍，也能给出准确的ANI估计值。

5. 数据库灵活：提供预构建的细菌、病毒、真核生物数据库，也支持用户轻松构建自定义数据库（如使用自己的MAGs）。

6. 读长兼容性好：同时适用于短读长和长读长测序数据，在牛津纳米孔的独立基准测试中表现最佳。

简而言之，Sylph 是一个旨在快速、精准地解析宏基因组样本物种组成和基因组相似性的强大工具。


##  安装sylph

#### 方式 1: conda install 
```sh
conda install -c bioconda sylph
```

#### 方式 2: 预编译可行性程序 (x86-64 linux statically compiled executable)

如果你使用的是x86-64系统，可以直接下载二进制文件并直接使用，无需安装。

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

注意：该二进制文件使用了不同的libraries库集（musl替代glibc）进行编译，这可能会影响性能。

## 快速开始

#### 进行宏基因组样本谱比对[GTDB-R226](https://gtdb.ecogenomic.org/)（包含143,614个代表性细菌/古菌物种基因组）

```sh
conda install -c bioconda sylph

# download GTDB-R226 pre-built database (~18.4 GB)
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb

# multi-sample paired-end profiling (sylph version >= 0.6)
sylph profile gtdb-r226-c200-dbv1.syldb -1 *_1.fastq.gz -2 *_2.fastq.gz -t 8 > profiling.tsv

# multi-sample single-end profiling
sylph profile gtdb-r226-c200-dbv1.syldb *.fastq -t 8 > profiling.tsv
```


## 教程、手册与预编译数据库

### 预编译数据库

预编译数据库可以从官网下载 [pre-built-databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/).

### [Cookbook](https://sylph-docs.github.io/sylph-cookbook/)

For common use cases and fast explanations, see the above [cookbook](https://sylph-docs.github.io/sylph-cookbook/).

### Tutorials
1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://sylph-docs.github.io/5%E2%80%90minute-sylph-tutorial/)
2. #### [Taxonomic profiling against GTDB database with MetaPhlAn-like output format](https://github.com/bluenote-1577/sylph/wiki/Taxonomic-profiling-with-the-GTDB%E2%80%90R214-database)

### Manuals
1. #### [Output format (TSV) and containment ANI explanation](https://github.com/bluenote-1577/sylph/wiki/Output-format)
2. #### [Taxonomic integration and custom taxonomies](https://github.com/bluenote-1577/sylph/wiki/Incorporating-taxonomic-information-into-sylph-with-sylph%E2%80%90tax)

### [sylph-tax](https://github.com/bluenote-1577/sylph-tax) 

To incorporate *taxonomy* into sylph's outputs, see the [sylph-tax repository](https://github.com/bluenote-1577/sylph-tax). 

> [!TIP] 
> The new [sylph-tax](https://github.com/bluenote-1577/sylph-tax) program replaces the old [sylph-utils](https://github.com/bluenote-1577/sylph-utils) repository. 

## Changelog

#### Version v0.8.0 - 2024-12-12. 

* Made the `inspect` option much less memory intensive. Slightly changed outputs when no genomes are found.

See the [CHANGELOG](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md) for complete details.

## Citing sylph

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.
