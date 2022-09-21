rm(list = ls())
library(WGCNA)

# 读取基因表达矩阵数据
tpm <- read.csv("批次效应output.csv", row.names=1)
head(tpm)

####数据清洗，筛除缺失值10%以上
library(DMwR) #用自带数据集algae，18个变量，200个观测值
library(VIM)

sum(!complete.cases(tpm)) #查看含有缺失值的样本个数
#直接删除所有含缺失值的样本
fpkm1<-na.omit(tpm) 
sum(!complete.cases(algae1))
#只删除缺失值过多的样本：缺失值个数大于列数的20%
algae2<-fpkm[-manyNAs(fpkm,0.1),] #数据框的“删除行”操作
sum(!complete.cases(algae2))
fpkm=algae2
### 选取基因方法 ###

## 第一种，通过标准差选择
## 计算每个基因的标准差
tpm_sd = apply(tpm,1,sd)#1是对每一行，2是对每一列
## 使用标准差对基因进行降序排序
tpm_sd_sorted = order(tpm_sd, decreasing = T)
## 选择前5000个标准差较大的基因
tpm_num = tpm_sd_sorted[1:5000]
## 从表达矩阵提取基因
tpm_filter = tpm[tpm_num,]
## 对表达矩阵进行转置
WGCNA_matrix = t(tpm_filter)#变成行名是样本，列名是基因
## 保存过滤后的数据
save(WGCNA_matrix, file = "data/Step01-tpm_sd_filter.Rda")

## 第二种，使用绝对中位差选择，推荐使用绝对中位差
WGCNA_matrix = t(tpm[order(apply(tpm,1,mad), decreasing = T)[1:5000],])#mad代表绝对中位差
save(WGCNA_matrix, file = "data/Step01-tpm_mad_filter-hp.Rda")

## 第三种，使用全基因
WGCNA_matrix = t(tpm)
save(WGCNA_matrix, file = "data/Step01-tpm_allgene.Rda")


### 去除缺失值较多的基因/样本 ###
rm(list = ls())
# 加载表达矩阵
load(file = "data/Step01-tpm_mad_filter-hp.Rda")
# 加载WGCNA包
library(WGCNA)
# 判断是否缺失较多的样本和基因
datExpr0 = WGCNA_matrix
gsg = goodSamplesGenes(datExpr0, verbose = 3)

# 是否需要过滤，TRUE表示不需要，FALSE表示需要
gsg$allOK
# 当gsg$allOK为TRUE，以下代码不会运行，为FALSE时，运行以下代码过滤样本
if (!gsg$allOK)
{
  # 打印移除的基因名和样本名
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 提取保留的基因和样本
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

### 通过样本聚类识别离群样本，去除离群样本 ###
sampleTree = hclust(dist(datExpr0), method = "average");#使用hclust函数进行均值聚类
# 绘制样本聚类图确定离群样本
sizeGrWindow(30,9)
pdf(file = "figures/Step01-sampleClustering.pdf", width = 30, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# 根据上图判断，需要截取的高度参数h
abline(h = 62, col = "red")#在120的地方画条线
dev.off()

# 去除离群得聚类样本，cutHeight参数要与上述得h参数值一致
clust = cutreeStatic(sampleTree, cutHeight = 62, minSize = 10)
table(clust)
# clust
# 0   1  0就是要去除的，1就是要保存的
# 15 162

# clust 1聚类中包含着我们想要得样本，将其提取出来
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr)#基因数
nSamples = nrow(datExpr)#样本数
save(datExpr, nGenes, nSamples,file = "data/Step01-WGCNA_input.Rda")

