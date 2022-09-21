# 导入数据及数据预处理

rows <- read.delim("~/Desktop/杨舒实验结果/A02-2022_07_09公共数据集生信/8.免疫细胞浸润要改/13.immcor/CIBERSORT-Results_filter.txt")
cols <- read.delim("~/Desktop/杨舒实验结果/A02-2022_07_09公共数据集生信/8.免疫细胞浸润要改/13.immcor/LASSOTURN.txt")
row.names(rows)<-rows[,1]#修改行名
rows<-rows[,-1]
row.names(cols)<-cols[,1]#修改行名
cols<-cols[,-1]
# 构建相关关系矩阵
library(psych)
data.corr <- corr.test(rows, cols, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值
# 画热图
library(pheatmap)
pheatmap(data.r, clustering_method="average", cluster_rows=F)
data.r.fmt <- matrix(sprintf("%.2f", data.r), nrow=nrow(data.p))  # 只保留小数点后两位
pheatmap(data.r, clustering_method="average", cluster_rows=F, display_numbers=data.r.fmt)
#加星号
getSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
  }
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)
pheatmap(data.r, clustering_method="average", cluster_rows=F, display_numbers=sig.mat)


