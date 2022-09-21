# 清空环境变量
rm(list = ls())
#加载R包
library(WGCNA)
# 允许多线程运行enableWGCNAThreads()
#mac不能多线程运行，只能在代码里加一句disableWGCNAThreads()
# 加载表达矩阵
load("data/Step01-WGCNA_input.Rda")

## 选择软阈值
sizeGrWindow(9,5);
par(mfrow = c(1,2));
powers = c(c(1:10), seq(from = 12, to=20, by=1))
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers)[[2]]
cex1=1.5
pdf(file="figures/Step03-softThresholding.pdf",width = 14)
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",main = paste("Scale independence"))
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], 
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")
plot(RpowerTable[,1], RpowerTable[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(RpowerTable[,1], 
     RpowerTable[,5], 
     labels=powers, 
     cex=cex1,col="red")
dev.off()

## 无尺度网络验证
softpower=5

ADJ = abs(cor(datExpr,use="p"))^softpower
k = as.vector(apply(ADJ,2,sum,na.rm=T))

pdf(file = "figures/Step03-scaleFree.pdf",width = 14)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k,main="Check scale free topology")
dev.off()

## 计算邻接矩阵
adjacency = adjacency(datExpr,power=softpower)
## 计算TOM拓扑矩阵
TOM = TOMsimilarity(adjacency)
## 计算相异度
dissTOM = 1- TOM

#模块初步聚类分析
library(flashClust)
geneTree = flashClust(as.dist(dissTOM),method="average")
#绘制层次聚类树
pdf(file = "figures/Step03-GeneClusterTOM-based.pdf")
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based",labels=FALSE,hang=0.04)
dev.off()

#构建初步基因模块
#设定基因模块中至少30个基因
minModuleSize=30
# 动态剪切树识别网络模块
dynamicMods = cutreeDynamic(dendro = geneTree,#hclust函数的聚类结果
                            distM = dissTOM,#
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)#设定基因模块中至少30个基因
# 将标签转换为颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)#看聚类到哪些模块，哪些颜色

pdf(file="figures/Step03-DynamicTreeCut.pdf",width=9,height=5)
plotDendroAndColors(dendro = geneTree, 
                    colors = dynamicColors, 
                    groupLabels = "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE,
                    main = "Gene dendrogram and module colors")
dev.off()


#前面用动态剪切树聚类了一些模块，现在要对这步结果进一步合并，合并相似度大于0.75的模块，降低网络的复杂度
#计算基因模块特征向量
MEList = moduleEigengenes(datExpr,colors=dynamicColors)#计算特征向量
MEs = MEList$eigengenes;#提取特征向量
MEDiss1 = 1-cor(MEs);#计算相异度
METree1 = flashClust(as.dist(MEDiss1),method="average")#对相异度进行flashClust聚类
#设置特征向量相关系数大于0.75
MEDissThres = 0.25;#相异度在0.25以下，也就是相似度大于0.75，对这些模块合并
#合并模块
merge = mergeCloseModules(datExpr, #合并相似度大于0.75的模块
                          dynamicColors,
                          cutHeight = MEDissThres, 
                          verbose=3)
mergedColors = merge$colors

table(dynamicColors)#动态剪切树的模块颜色
table(mergedColors)#合并后的模块颜色，可以看到从18个模块变成了14个模块

mergedMEs = merge$newMEs#合并后的14个模块

#重新命名合并后的模块
moduleColors = mergedColors;
colorOrder = c("grey",standardColors(50));
moduleLabels = match(moduleColors,colorOrder)-1;
MEs = mergedMEs;
MEDiss2 = 1-cor(MEs);#计算相异度
METree2 = flashClust(as.dist(MEDiss2),method="average");#对合并后的模块进行聚类

#绘制聚类结果图
pdf(file="figures/Step03-MECombined.pdf",width=12,height=5)
par(mfrow=c(1,2))
plot(METree1,xlab="",sub="",main="Clustering of ME before combined")# METree1是动态剪切树形成的模块
abline(h=MEDissThres,col="red")#相异度为0.25
plot(METree2,xlab="",sub="",main="Clustering of ME after combined")# METree2是合并后的模块
dev.off()

pdf(file="figures/Step03-MergedDynamics.pdf",width=10,height=4)
plotDendroAndColors(dendro = geneTree,#剪切树
                    colors = cbind(dynamicColors,mergedColors),#将两种方法形成的模块颜色合并在一起
                    groupLabels = c("Dynamic Tree Cut","Merged Dynamics"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    main="Gene Dendrogram and module colors")
dev.off()

# 模块中基因数
write.table(table(moduleColors),"data/Step03-MEgeneCount.txt",quote = F,row.names = F)

# 保存构建的网络信息
moduleColors=mergedColors
colorOrder=c("grey", standardColors(50))
moduleLabels=match(moduleColors, colorOrder)-1
MEs=mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file="data/Step03-Step_by_step_buildnetwork.rda")




## 绘制样本间的相关性
load("data/Step01-WGCNA_input.Rda")
load(file="data/Step03-Step_by_step_buildnetwork.rda")


MEs = orderMEs(MEs)

sizeGrWindow(5,7.5);
pdf(file = "figures/Step03-moduleCor.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

## TOMplot

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softpower); 

nSelect = 400
# 随机选取400个基因进行可视化，s设置seed值，保证结果的可重复性
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# 对选取的基因进行重新聚类
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = mergedColors[select];
# 打开一个绘图窗口
sizeGrWindow(9,9)
pdf(file = "figures/Step03-TOMplot.pdf", width = 9, height = 9);
# 美化图形的设置
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes")
dev.off()

