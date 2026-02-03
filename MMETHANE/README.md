# MMETHANE：微生物和代谢物宿主分析引擎

MMETHANE（Microbes and METabolites to Host Analysis Engine）是一个计算工具，能够从联合代谢组学和微生物测序数据中做出可解释的二元结果预测。它支持 16S rRNA 扩增子测序或宏基因组学数据，输出英语语言解释其决策过程，并提供用于决策的基础数据可视化。

**主要特性：**
- 支持多种机器学习模型：MMETHANE（自定义模型）、前馈神经网络、LASSO 逻辑回归、随机森林、AdaBoost
- 可处理 16S rRNA 扩增子测序和 WGS/宏基因组学数据
- 集成代谢物数据分析和微生物序列数据分析
- 提供可解释的预测结果和可视化
- 支持数据预处理和特征选择
- 多线程并行处理支持

## 安装

建议在 Python 3.11 的虚拟环境中安装 MMETHANE 及其所有依赖项。

### 创建虚拟环境
```bash
conda create -n mmenv python=3.11
conda activate mmenv
```

### 安装 MMETHANE

**方法1：使用 pip 安装**
```bash
pip install mmethane==0.2
```

**方法2：从源码安装**
```bash
git clone https://github.com/gerberlab/mmethane.git
cd mmethane
pip install -r requirements.txt
pip install -r requirements_c.txt
conda install -c conda-forge ete3 -y
```

## 快速开始

要尝试 MMETHANE，可以使用提供的示例配置文件：

### 1. 运行 MMETHANE 模型
```bash
mmethane -c config_files/sample.cfg -o <绝对路径到输出文件夹/>
```

### 2. 运行 LASSO 逻辑回归模型
```bash
mmethane -c config_files/sample_LR.cfg -o <绝对路径到输出文件夹/>
```

运行这些命令将：
1. 处理示例数据集
2. 运行指定的模型
3. 输出包含可视化的 HTML 文件

可视化文件将位于：`<输出文件夹>/mmethane_franzosa/seed_0/visualization.html`

## 配置文件说明

MMETHANE 使用配置文件来指定所有输入参数和运行选项。配置文件采用 INI 格式，包含多个部分：

### [description] 部分
- `tag`：用于数据集文件夹的名称
- `in_path`：（可选）包含所有输入数据的文件夹路径
- `process_data`：（可选）是否处理数据，默认为 True
- `run_model`：（可选）是否运行模型，默认为 True

### [run] 部分
- `model`：**必需**，要运行的模型，选项："mmethane"、"ffn"、"lr"、"rf"、"adaboost"
- `seed`：**必需**，运行的随机种子
- `run_name`：（可选）输出文件夹中的日志文件夹名称，如未指定则使用 "tag"
- `out_path`：（可选）日志输出路径
- `dtype`：（可选）数据类型，如 "metabs, otus"
- `parallel`：（可选）并行线程数

### [data] 部分
- `subject_data`：**必需**，受试者元数据文件路径（CSV/TSV），行是样本，列至少包含结果变量
- `outcome_variable`：**必需**，受试者数据中结果变量的列名
- `outcome_positive_value`：（可选）如果结果值不是整数（0 和 1），指定表示阳性结果的字符串
- `outcome_negative_value`：（可选）指定阴性结果值，或与 `outcome_positive_value` 一起使用以排除具有第三种结果的样本
- `covariate_variable`：（可选）指定数据应拆分依据的协变量变量
- `sample_id_column`：（可选）如果样本 ID 不是受试者数据的索引，指定样本 ID 列

### [sequence_data] 部分
- `data_type`：**必需**，数据类型，选项："16s" 或 "WGS"
- `data`：**必需**，序列数据文件路径，格式为样本行和序列列
- `reference_tree`：WGS 数据必需，用于处理数据的参考树
- `sequences`：如果序列数据列标签不是实际的 RNA 序列字符串，则需要此文件
- `taxonomy`：（可选）输入分类文件
- `tree`：（可选）数据中序列的系统发育树
- `distance_matrix`：（可选）距离矩阵
- `samples_dimension`：（可选）样本维度，选项："rows" 或 "columns"，默认为 "rows"

### [sequence_preprocessing] 部分（可选）
- `process_before_training`：（可选）是否在训练前处理数据，默认为 False
- `percent_present_in`：（可选）微生物必须高于检测限的样本百分比
- `limit_of_detection`：（可选）检测限，默认为 0
- `cov_percentile`：（可选）变异系数百分位数
- `transformations`：（可选）转换，如 "relative_abundance"

### [metabolite_data] 部分
- `data`：**必需**，代谢物数据文件路径
- `meta_data`：必需（除非使用原始 CDI 数据），包含 HMDB、KEGG 或 InChIKey 标识符的元数据
- `taxonomy`：（可选）代谢物的 ClassyFire 分类文件
- `fingerprint_type`：（可选）指纹类型，默认为 "pubchem"
- `similarity_matrix`：（可选）相似性矩阵路径
- `distance_matrix`：（可选）距离矩阵路径
- `samples_dimension`：（可选）样本维度，默认为 "rows"

### [metabolite_preprocessing] 部分（可选）
- `process_before_training`：（可选）是否在训练前处理数据，默认为 False
- `percent_present_in`：（可选）代谢物必须高于检测限的样本百分比
- `limit_of_detection`：（可选）检测限，默认为 0
- `cov_percentile`：（可选）变异系数百分位数
- `transformations`：（可选）转换，如 "log,standardize"

## 示例配置文件

### 基本 MMETHANE 配置示例
```ini
[description]
tag:mmethane_franzosa
in_path:datasets/FRANZOSA/
process_data:True
run_model:True

[run]
model:MMETHANE
seed: 0
dtype: metabs, otus
parallel: 6

[data]
subject_data:${description:in_path}/metadata_cv.csv
outcome_variable:Study.Group
sample_id_column:Sample
outcome_negative_value:CD,UC
outcome_positive_value:Control

[sequence_data]
samples_dimension:columns
data_type:WGS
data: ${description:in_path}/merged_sp_cts.csv
reference_tree:mmethane/utilities/phylo_references/mpa_v31_CHOCOPhlAn_201901_species_tree.nwk.txt

[metabolite_data]
samples_dimension:rows
data:${description:in_path}/mtb.tsv
meta_data:${description:in_path}/mtp_map_wInchiKey.csv
taxonomy:${description:in_path}/classy_fire_df.csv
fingerprint_type:pubchem
skip_taxonomy:True

[metabolite_preprocessing]
transformations:log,standardize
percent_present_in:15
limit_of_detection:0
```

### LASSO 逻辑回归配置示例
```ini
[description]
tag:mmethane_franzosa
in_path:../datasets/FRANZOSA/
process_data:False
run_model:True

[run]
model:lr
out_path: logs/
run_name:lr_franzosa
seed: 0
dtype: metabs, otus
parallel: 6
use_ray: 0

[data]
subject_data:${description:in_path}/metadata_cv.csv
outcome_variable:Study.Group
sample_id_column:Sample
outcome_negative_value:CD,UC
outcome_positive_value:Control

[sequence_data]
samples_dimension:columns
data_type:WGS
data: ${description:in_path}/merged_sp_cts.csv
reference_tree:./utilities/phylo_references/mpa_v31_CHOCOPhlAn_201901_species_tree.nwk.txt
```

## 输入数据格式

### 序列数据
- **16S 数据**：ASV/OTU 计数表，行为样本，列为序列特征
- **WGS 数据**：物种丰度表，支持多种格式
- 序列可以用标识符或实际的 RNA 序列字符串命名

### 代谢物数据
- 代谢物丰度表，行为样本，列为代谢物特征
- 需要元数据映射文件，包含 HMDB、KEGG 或 InChIKey 标识符
- 支持 PubChem、RDKit、Morgan 和 MQN 指纹

### 受试者元数据
- CSV/TSV 格式，包含样本结果和可能的协变量信息
- 样本作为行，结果变量作为列
- 支持二元和多类结果

## 支持的模型

1. **MMETHANE**：自定义神经网络模型，专门为多组学数据设计
2. **FFN（前馈神经网络）**：标准神经网络模型
3. **LR（LASSO 逻辑回归）**：带 L1 正则化的逻辑回归
4. **RF（随机森林）**：集成树模型
5. **AdaBoost（自适应梯度提升）**：提升算法

## 输出和可视化

MMETHANE 生成以下输出：
- **HTML 可视化文件**：包含模型决策的可视化解释
- **日志文件**：详细的运行日志和性能指标
- **处理后的数据**：预处理后的数据集
- **模型检查点**：训练好的模型权重

可视化包括：
- 特征重要性分析
- 数据分布图
- 模型性能指标
- 预测结果解释

## 依赖项

主要依赖包括：
- PyTorch（2.5.1）和 PyTorch Lightning（2.4.0）
- DGL（2.2.0）和 DGLLife（0.3.2）用于图神经网络
- scikit-learn（1.5.2）用于传统机器学习模型
- RDKit（2024.3.6）用于化学信息学
- ete3（3.1.3）用于系统发育树处理
- PubChemPy（1.0.4）用于代谢物数据访问
- 其他标准科学计算库（numpy, pandas, scipy, matplotlib, seaborn）

## 许可证

GPL-3.0 许可证

## 参考文献

相关出版物即将发布。请关注 GitHub 仓库获取最新信息。

**GitHub 仓库**：https://github.com/gerberlab/mmethane

**创建时间**：2024年11月8日

**最后更新**：2025年9月19日

**主要语言**：Jupyter Notebook、Python

**作者**：Jennifer Dawkins（Gerber Lab）