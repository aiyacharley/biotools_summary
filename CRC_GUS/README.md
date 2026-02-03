# CRC_GUS：结直肠癌肠道微生物 β-葡萄糖醛酸酶图谱

CRC_GUS 是一个综合性研究工具包，用于构建结直肠癌（CRC）患者肠道微生物 β-葡萄糖醛酸酶（gmGUSs）的图谱，并分析微生物-gmGUS-代谢物轴在结直肠癌发生发展中的作用。该项目包含从 GUS 基因鉴定到完整数据分析的完整工作流程。

**主要特性：**
- 完整的 GUS 鉴定流程：基于相似性、结构域和保守残基
- Loop 分类系统：识别 Mini-Loop2 等关键结构特征
- 分类学注释数据库构建：支持 Kraken2 分类
- 全面的 R 分析流程：从基础统计到机器学习模型
- 多队列验证：在独立队列中验证发现
- 代谢物关联分析：揭示功能生物学联系

## 项目概述

该研究构建了来自公共 CRC 队列的 550 个 gmGUSs 的图谱，使用了 114 个参考 GUS、三个 GUS 结构域和七个保守残基。研究发现：
- 晚期 CRC 中 Mini-Loop2 和 GUS 携带物种（Bacteroides cellulosilyticus 和 Bacteroides nordii）富集
- 38 个差异 gmGUSs 有效区分患者与对照组（AUC > 0.8）
- 基于 B. cellulosilyticus 的五个 gmGUSs 的 GUSscore 模型能很好预测 CRC 结局
- 微生物-gmGUS-代谢物轴在氨基酸代谢、维生素生物合成、细菌行为和 LPS 生物合成中的特异性生物学联系

## 项目结构

```
CRC_GUS/
├── 01.GUS_Identification/     # GUS 鉴定流程
│   ├── testInput/            # 测试输入文件
│   ├── resource/             # 参考数据库和结构域
│   ├── script/               # 核心 Perl 脚本
│   ├── work.sh              # 工作流程脚本
│   └── README.md
├── 02.Loop_Classification/   # Loop 分类流程
│   ├── testInput/
│   ├── resource/
│   ├── script/
│   ├── work.sh
│   └── README.md
├── 03.Taxonomic_Annotation/  # 分类学注释数据库构建
│   ├── krakenDB.pl
│   └── README.md
├── 04.DataAnalysis/         # 数据分析 R 脚本
│   ├── 00.rawdata/          # 原始数据文件
│   ├── 01.rarefaction_barplot.R
│   ├── 02.LoopAnalysis.R
│   ├── 03.DiversityAnalysis.R
│   ├── 04.speciesAnalysis.R
│   ├── 05.GUSanalysis.R
│   ├── 06.RFmodel.R
│   ├── 07.GUSscoreModel.R
│   ├── 08.speciesCorr.R
│   ├── 09.metaboliteAnalysis.R
│   ├── 10.EXPandRNAseq.R
│   └── README.md
├── LICENSE
└── README.md
```

## 系统要求

### GUS 鉴定和分类流程
- **Perl 5**：版本 26.2（v5.26.2）
- **BLAST+**：blastp 2.12.0+
- **HMMER**：hmmsearch 3.3.2
- **Clustal Omega**：1.2.3（用于 loop 分类）
- **Kraken2**：2.0.7-beta（用于分类学注释）
- **操作系统**：Linux（macOS 理论上可行）

### 数据分析流程
- **R**：版本 4.3.2 或更高
- **R 包**：vegan, dplyr, ggplot2, ggpubr, viridis, reshape2, ape, VennDiagram, ggtree, treeio, pheatmap, psych, coin, caret, pROC, Boruta, randomForest, survival, survminer, glmnet, timeROC, readxl, clusterProfiler, org.Hs.eg.db, grid, ggsci, ggsignif, ggbreak, Rmisc, ggprism, EnhancedVolcano, tibble
- **操作系统**：macOS（Linux 理论上可行）

## 安装

### 1. 克隆仓库
```bash
git clone https://github.com/yr2008-UM/CRC_GUS.git
cd CRC_GUS
```

### 2. 安装 Perl 依赖（GUS 鉴定流程）
```bash
# 确保已安装所需软件
# Ubuntu/Debian
sudo apt-get install perl ncbi-blast+ hmmer clustalo kraken2

# 或从源码安装
# BLAST: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
# HMMER: http://hmmer.org/download.html
# ClustalO: http://www.clustal.org/omega/
# Kraken2: https://ccb.jhu.edu/software/kraken2/index.shtml
```

### 3. 安装 R 依赖（数据分析流程）
```r
# 安装所需 R 包
install.packages(c("vegan", "dplyr", "ggplot2", "ggpubr", "viridis", "reshape2",
                   "ape", "VennDiagram", "ggtree", "treeio", "pheatmap", "psych",
                   "coin", "caret", "pROC", "Boruta", "randomForest", "survival",
                   "survminer", "glmnet", "timeROC", "readxl", "clusterProfiler",
                   "org.Hs.eg.db", "grid", "ggsci", "ggsignif", "ggbreak", "Rmisc",
                   "ggprism", "EnhancedVolcano", "tibble"))

# 如果失败，尝试使用 BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggtree", "treeio", "clusterProfiler", "org.Hs.eg.db"))
```

## 使用示例

### 1. GUS 鉴定流程

#### 运行完整工作流程
```bash
cd 01.GUS_Identification
bash work.sh
```

#### 逐步执行
```bash
# 步骤1：对单基因进行 blastp 和 hmmsearch
blastp -query testInput/InputGeneCatalogue.prot.fa -db resource/references \
  -evalue 0.05 -out 01.blastp.result.xls -outfmt 6 -num_threads 4

hmmsearch --tblout 02.hmm.result.tblout --domtblout 02.hmm.result.domtblout \
  --pfamtblout 02.hmm.result.pfamtblout -o 02.hmm.result.xls -E 0.05 --cpu 4 \
  resource/domains.hmm testInput/InputGeneCatalogue.prot.fa

# 步骤2：基于 blastp 和 hmmsearch 结果筛选
perl -ne 'my@or=split/\s+/;print if $or[2] > 25 & $or[10] < 0.05 & !$check{$or[0]};$check{$or[0]}=1;' \
  01.blastp.result.xls > 01.blastp.result.screen.xls

perl -ne 'BEGIN{$/="\n>";for my $l  (`less testInput/InputGeneCatalogue.prot.fa`){chomp$l;$l=~s/^>//;my $id=$1 if $l=~/^(\S+)/;$seq{$id}=$l;}$/="\n";for my $l (`less 01.blastp.result.screen.xls`){my@or=split/\s+/,$l;$check{$or[0]}=1;}}next if $_=~/^#/;my@or=split/\s+/;$hash{$or[0]}{$or[2]}=1 if $or[4] < 0.05;END{foreach my $id (sort keys %hash){print ">$seq{$id}\n" if $hash{$id}{"Glyco_hydro_2_C"} && $hash{$id}{"Glyco_hydro_2_N"} && $hash{$id}{"Glyco_hydro_2"} && $check{$id} }}' \
  02.hmm.result.tblout > 03.1.hmm.and.blast.screen.fa

# 步骤3：基于 motif 进一步筛选候选序列
perl script/01.checkMotifs.pl 03.1.hmm.and.blast.screen.fa 03.2.checkMotifs.hmm.fa
```

### 2. Loop 分类流程
```bash
cd 02.Loop_Classification

# 多序列比对
clustalo -i testInput/input.fa -o 04.clustalO.output.txt \
  --threads 4 --p1 resource/12refsForLoop.clustalO.aln.txt

# Loop 分类
perl script/02.loopClass.pl 04.clustalO.output.txt 05.loop.stat.xls > 06.loopclass.log
```

### 3. 分类学注释数据库构建
```bash
cd 03.Taxonomic_Annotation

# 下载必要文件（需提前下载）
# nt.gz: https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
# nucl_gb.accession2taxid: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# taxdump.tar.gz: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# 构建标准数据库
kraken2-build --standard --threads 5 --db dirForDB

# 提取细菌、古菌和病毒的核苷酸序列
perl krakenDB.pl nodes.dmp names.dmp nucl_gb.accession2taxid nt.fna dirForDB/nt/library.fna

# 构建基因组文库
kraken2-build --download-library nt --db dirForDB --threads 10
kraken2-build --build --db dirForDB --threads 10
```

### 4. 数据分析流程（R 脚本）

#### 设置工作目录和加载数据
```r
# 设置工作目录
setwd("path/to/CRC_GUS/04.DataAnalysis")

# 加载所需库
library(vegan)    # 多样性分析
library(dplyr)    # 数据操作
library(ggplot2)  # 可视化

# 读取数据
group <- read.csv("00.rawdata/group.csv", header = TRUE, row.names = 1)
GUSabun_ABS <- read.csv("00.rawdata/GUSabun_ABS.csv", header = TRUE, row.names = 1)
referenceStat <- read.csv("00.rawdata/referenceStat.csv", header = TRUE, row.names = 1)
GUSsStat <- read.csv("00.rawdata/GUSsStat.csv", header = TRUE, row.names = 1)
```

#### 稀释曲线分析（01.rarefaction_barplot.R）
```r
# 生成所有样本的稀释曲线
data.spec <- specaccum(t(GUSabun_ABS), method = "random")
plot(data.spec, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0,
     ci.col = "lightblue", xlab = "样本数", ylab = "gmGUS 数量")

# 添加健康组曲线
healthy_samples <- group %>% dplyr::filter(Stage == "Healthy") %>% rownames()
data.spec <- specaccum(t(GUSabun_ABS[, healthy_samples]), method = "random")
plot(data.spec, ci.type = "poly", col = '#cecccb', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", add = TRUE)

# 添加 CRC 组曲线
crc_samples <- group %>% dplyr::filter(Stage %in% c('S0','SI_II','SIII_IV')) %>% rownames()
data.spec <- specaccum(t(GUSabun_ABS[, crc_samples]), method = "random")
plot(data.spec, ci.type = "poly", col = '#f3764a', lwd = 2, ci.lty = 0,
     ci.col = "lightblue", add = TRUE)
```

#### GUS 差异分析（05.GUSanalysis.R）
```r
# 组间统计函数
group_mean <- function(group.subgroup, data, colInfo = "Group"){
  tmpg <- unique(group.subgroup[, colInfo])
  tmpn <- length(tmpg)
  rows <- nrow(data)

  result <- matrix(nrow = rows, ncol = tmpn * 3)
  name <- NULL

  for(i in tmpg){
    name <- c(name, paste(i, 'Median', sep = '_'),
                    paste(i, 'Mean', sep = '_'),
                    paste(i, 'SD', sep = '_'))
  }
  colnames(result) <- name

  for(i in tmpg){
    tmpdata <- data[, rownames(group.subgroup[which(group.subgroup[, colInfo] == i), , drop = F])]
    name <- c(paste(i, 'Median', sep = '_'), paste(i, 'Mean', sep = '_'), paste(i, 'SD', sep = '_'))

    for (j in 1:rows){
      tmp.median <- median(as.numeric(tmpdata[j, ]))
      tmp.mean <- mean(as.numeric(tmpdata[j, ]))
      tmp.SD <- sd(as.numeric(tmpdata[j, ]))
      result[j, name] <- c(tmp.median, tmp.mean, tmp.SD)
    }
  }

  rownames(result) <- rownames(data)
  result <- data.frame(result)
  return(result)
}

# Wilcoxon 检验和 FDR 校正
wilcoxon.FDR.TAX <- function(data, group){
  data <- data[, rownames(group)]
  tmpg <- unique(group$Stage)
  species.abun <- data

  result <- group_mean(group, data, 'Stage')

  # wilcoxon 检验
  for (i in c('Healthy')) {
    for (j in c('MP', 'S0', "SI_II", "SIII_IV")){
      sample1 <- rownames(group[which(group[, 'Stage'] == i), , drop = F])
      sample2 <- rownames(group[which(group[, 'Stage'] == j), , drop = F])

      wilcox.01.p <- apply(data, 1, function(x) {
        if(sum(x[c(sample1, sample2)]) == 0){
          return('NA')
        } else {
          test <- wilcox.test(x[sample1], x[sample2], conf.int = TRUE)
          p <- test$p.value
          return(p)
        }
      }) %>% as.numeric()

      wilcox.01.q <- rep("NA", length(wilcox.01.p))
      wilcox.01.q[which(wilcox.01.p < 0.05)] <- p.adjust(wilcox.01.p[which(wilcox.01.p < 0.05)], method = "fdr")

      wilname.p <- data.frame(p = wilcox.01.p, q = wilcox.01.q)
      colnames(wilname.p) <- paste0(c("wilcox.test.p", "wilcox.test.q"), paste0("(", i, ' VS ', j, ")"))
      result <- cbind(result, wilname.p)
    }
  }

  return(result)
}
```

#### 随机森林分类器（06.RFmodel.R）
```r
# 加载所需库
library(caret)      # 分类和回归训练
library(pROC)       # ROC 分析
library(Boruta)     # 特征选择
library(randomForest) # 随机森林

# 特征选择使用 Boruta 算法
boruta_output <- Boruta(Group ~ ., data = training_data, doTrace = 2)
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)

# 模型训练
rf_model <- randomForest(
  Group ~ .,
  data = training_data[, c("Group", boruta_signif)],
  ntree = 500,
  mtry = sqrt(length(boruta_signif)),
  importance = TRUE
)

# 模型评估
predictions <- predict(rf_model, testing_data)
confusion_matrix <- confusionMatrix(predictions, testing_data$Group)

# ROC 分析
roc_curve <- roc(testing_data$Group, as.numeric(predictions))
auc_value <- auc(roc_curve)

# 变量重要性
var_importance <- importance(rf_model)
varImpPlot(rf_model)
```

#### GUSscore 模型（07.GUSscoreModel.R）
```r
# 加载生存分析库
library(survival)    # 生存分析
library(survminer)   # 生存可视化
library(glmnet)      # LASSO 回归
library(timeROC)     # 时间依赖的 ROC

# Cox 比例风险回归
cox_model <- coxph(Surv(time, status) ~ ., data = survival_data)

# LASSO 回归选择特征
x <- as.matrix(survival_data[, -c(1:2)])  # 排除时间和状态列
y <- Surv(survival_data$time, survival_data$status)

cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1)
best_lambda <- cv_fit$lambda.min

lasso_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)
selected_features <- coef(lasso_model)@Dimnames[[1]][which(coef(lasso_model) != 0)]

# 构建 GUSscore
GUSscore <- apply(x[, selected_features], 1, function(row) sum(row * coef(lasso_model)[selected_features]))

# 生存分析
surv_fit <- survfit(Surv(time, status) ~ ifelse(GUSscore > median(GUSscore), "High", "Low"),
                    data = survival_data)

ggsurvplot(surv_fit, data = survival_data, risk.table = TRUE,
           pval = TRUE, conf.int = TRUE, palette = c("#E7B800", "#2E9FDF"),
           legend.labs = c("High GUSscore", "Low GUSscore"))
```

## 输入数据格式

### 1. GUS 鉴定输入
- **蛋白质序列文件**：FASTA 格式的单基因蛋白质序列
- **参考数据库**：BLAST 格式的参考 GUS 序列
- **结构域 HMM 文件**：HMMER 格式的 GUS 结构域模型

### 2. 数据分析输入
- **样本元数据**（group.csv）：CSV 文件，行为样本 ID，列至少包含 Stage（Healthy, MP, S0, SI_II, SIII_IV）
- **GUS 丰度矩阵**（GUSabun_TPM.csv）：CSV 文件，行为 gmGUS ID，列为样本 ID，值为 TPM 归一化丰度
- **GUS 统计信息**（GUSsStat.csv）：CSV 文件，包含 gmGUS 的 loop 类别和分类学注释
- **参考序列统计**（referenceStat.csv）：CSV 文件，包含参考 GUS 的 loop 类别和分类学信息
- **系统发育树**（38gmGUS.tre）：Newick 格式的 38 个显著 gmGUSs 的系统发育树

## 输出结果

### 1. GUS 鉴定输出
- 筛选后的候选 GUS 序列（FASTA 格式）
- 长度统计文件
- 序列列表文件

### 2. 数据分析输出
- **图形输出**：所有论文中的图形（PDF/PNG 格式）
- **统计结果**：差异分析、相关性分析、模型性能指标
- **模型文件**：训练好的随机森林和 LASSO 模型
- **中间数据**：处理后的数据文件

## 数据分析流程概述

### 脚本执行顺序
1. **01.rarefaction_barplot.R**：稀释曲线和条形图（图 1b，补充图 1-2）
2. **02.LoopAnalysis.R**：loop 类别分布、丰度差异和分类组成（图 2a-c，补充图 3）
3. **03.DiversityAnalysis.R**：β-多样性（PCoA）和 α-多样性比较（图 2d-e）
4. **04.speciesAnalysis.R**：物种水平累计 GUS 丰度/数量和拷贝数变异（CNV）（图 2f-h，补充图 4）
5. **05.GUSanalysis.R**：gmGUSs 在 CRC 阶段的差异丰度分析和独立队列验证（图 3a，补充图 5, 7b）
6. **06.RFmodel.R**：CRC/腺瘤分类的随机森林分类器，包含特征选择和验证（图 3b-c，补充图 6-8）
7. **07.GUSscoreModel.R**：Cox 回归分析和 LASSO 构建的 GUSscore 模型预测 CRC 生存结局（图 3d-g，补充图 9-10）
8. **08.speciesCorr.R**：gmGUSs 与细菌物种的相关性分析（图 4，补充图 11）
9. **09.metaboliteAnalysis.R**：gmGUSs 与代谢物/KEGG Orthology 术语的相关性分析（图 5）
10. **10.EXPandRNAseq.R**：实验数据分析，包括酶测定、细胞实验和 RNA-seq 数据（图 6，补充图 12-14）

## 注意事项

1. **数据路径**：确保所有输入文件位于正确的路径（00.rawdata/ 目录）
2. **内存要求**：大规模数据分析可能需要足够的内存（建议 16GB+）
3. **并行处理**：部分脚本支持多线程，可调整线程数以提高性能
4. **外部依赖**：确保所有外部软件（BLAST、HMMER、Kraken2）正确安装并添加到 PATH
5. **R 包版本**：使用指定版本的 R 包以确保结果可重现

## 许可证

GNU General Public License v2.0 (GPL-2.0)

## 引用

使用此工具时，请引用相关的手稿（即将发布）。

**GitHub 仓库**：https://github.com/yr2008-UM/CRC_GUS

**创建时间**：2025年4月24日

**最后更新**：2025年9月9日

**主要语言**：R、Perl、Shell

**作者**：Junru Chen (junru.chen2019@hotmail.com) 等