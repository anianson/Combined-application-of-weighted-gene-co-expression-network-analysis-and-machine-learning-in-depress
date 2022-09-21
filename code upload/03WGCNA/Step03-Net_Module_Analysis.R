# 清空环境变量
rm(list = ls())
# 加载包
library(WGCNA)
# 加载表达矩阵
load("data/Step01-WGCNA_input.Rda")

# 读入临床信息
clinical = read.csv("sample_input.csv", row.names=1)
# 查看临床信息
head(clinical)
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)

pdf(file = "figures/Step04-Sample_dendrogram_and_trait_heatmap.pdf",
    width = 24);
# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    cex.dendroLabels = 1.5,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### 网络的分析
###### 基因模块与临床信息的关系
# 加载构建的网络
load(file = "data/Step03-Step_by_step_buildnetwork.rda")
# 对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")
write.table(file="data/Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="data/Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化。
pdf(file="figures/Step04-Module_trait_relationships.pdf",width=9,height=7)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=T,
               colors=blueWhiteRed(50),
               xColorOffset = strheight("M")/300, 
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=1,
               cex.lab=1,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()


##### 单一模块与某一表型相关性--GS_MM散点图
group = as.data.frame(datTraits$group)#改
# 分析自己感兴趣的临床信息，此处以M_stage为示例
names(group) = "group"#改
# 模块对应的颜色
modNames = substring(names(MEs), 3)

# 计算基因模块特征
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples));

# 对结果进行命名
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# 计算M分期基因特征显著性
geneTraitSignificance = as.data.frame(cor(datExpr,datTraits$group, 
                                          use = "p"));#改
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples));

# 对结果进行命名
names(geneTraitSignificance) = paste("GS.", names(group), sep="")#改
names(GSPvalue) = paste("p.GS.", names(group), sep="")#改

# 设置需要分析的模块名称，此处为brown模块
module = "pink"#改
# 提取brown模块数据
column = match(module, modNames);
moduleGenes = moduleColors==module;

# 可视化brown模块与M分期的相关性分析结果
sizeGrWindow(7, 7);
pdf(file="figures/Step04-Module_membership_vs_gene_significance-Immobility_time.pdf")
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Immobility_time",#改
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()


#即通过MM和GS进行筛选，这里的标准不固定，根据自己的数据进行更改。这里我用的标准是|GS|>0.2,|MM|>0.8
# 选择lightyellow模块
hub<- abs(geneModuleMembership$MMyellow)>0.8 & abs(geneTraitSignificance)>0.2
table(hub)
hub<-as.data.frame(dimnames(data.frame(datExpr))[[2]][hub]) 
write.csv(hub, "hubgene_GSMM_yellow_Immobility_time1.csv")#改

##########################GS_MM散点图绘制
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes, 1]), 
                   xlab = paste("Module Membership in", module, "module"), 
                   ylab = paste("Gene/significance for Immobility_time"),
                   main = paste("Module membership vs. gene significance"), 
                   pch = 20,
                   col="grey")
abline(h=0.2,v=0.8,col="red",lwd=1.5)
 module = "black"
 column = match(module, modNames)
 moduleGenes = moduleColors==module
 table(moduleGenes)
#moduleGenes
#FALSE  TRUE 
#38899  2726 
 black_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes]) 
 names(black_module)="genename"
  hub<- abs(geneModuleMembership$MMblack)>0.8 & abs(geneTraitSignificance)>0.2
  table(hub)
 #hub
 #FALSE  TRUE 
 #41563    62 
  hubgene_black <- dimnames(data.frame(datExpr))[[2]][hub]
   MM<-abs(geneModuleMembership[moduleGenes,column])
   GS<-abs(geneTraitSignificance[moduleGenes, 1])
   c<-as.data.frame(cbind(MM,GS))
   rownames(c)=black_module$genename
   head(c)
   # hub基因和module内全部基因进行匹配，匹配成功返回1，没有匹配到的返回0
    match<- black_module$genename %in% hubgene_black
   # 将匹配信息添加到散点图矩阵最后一列
    c$group<-match
    head(c)
     library(ggplot2)
     pdf("MM vs. GS_black_Immobility_time.pdf",width = 7,height = 7)
     ggplot(data=c, aes(x=MM, y=GS, color=group))+
       geom_point(size=1.5)+
       scale_colour_manual(values=c("grey60", "black"))+ 
       theme_bw()+
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+  
       labs(x="Module Membership in black module", 
            y="Gene significance for Immobility_time",
            title = "Module membership vs. gene significance,
        cor=0.38, p=3.2e−13 ")+
       theme(axis.title.x =element_text(size=14), 
             axis.title.y=element_text(size=14),
             axis.text = element_text(size = 12),
             axis.text.x = element_text(colour = "black"),
             axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
             plot.margin = unit(rep(2,4),'lines'))+
       theme(legend.position = 'none')+
       geom_hline(aes(yintercept=0.2),
                  colour="#5B9BD5",lwd=1,linetype=5)+
       geom_vline(aes(xintercept=0.8),colour="#5B9BD5",lwd=1,linetype=5)
     dev.off()
   
     
## 单一特征与所有模块关联分析
GS = as.numeric(cor(datTraits$group,datExpr, use="p"))#改
GeneSignificance = abs(GS);
pdf(file="figures/Step04-Gene_significance_for_group_across_module.pdf",width=9,height=5)#改
plotModuleSignificance(GeneSignificance, moduleColors, ylim=c(0,0.6), 
                       main="Gene significance for group across module");#改
dev.off()

## 模块中的hub基因
#### 为每一个模块寻找hub基因
HubGenes <- chooseTopHubInEachModule(datExpr,#WGCNA分析输入的表达矩阵
                                     moduleColors)#模块颜色信息
# 保存hub基因结果
write.table (HubGenes,file = "data/Step04-HubGenes_of_each_module.xls",quote=F,sep='\t',col.names = F)

#### 与某种特征相关的hub基因
NS = networkScreening(datTraits$group,#M分期
                      MEs,#
                      datExpr)#WGCNA分析输入的表达矩阵
# 将结果写入到文件
write.table(NS,file="data/Step04-Genes_for_group.xls",quote=F,sep='\t')

## 模块GO/KEGG分析
# 加载R包
library(anRichment)
library(clusterProfiler)
##### GO分析
# 构建GO背景基因集
GOcollection = buildGOcollection(organism = "rat")#改
geneNames = colnames(datExpr)
# 将基因SYMBOL转换为ENTREZID基因名
geneID = bitr(geneNames,fromType = "SYMBOL", toType = "ENTREZID", 
                 OrgDb = "org.Rn.eg.db", drop = FALSE)
# 将基因名对应结果写入文件中
write.table(geneID, file = "data/Step04-geneID_map.xls", sep = "\t", quote = TRUE, row.names = FALSE)
# 进行GO富集分析
GOenr = enrichmentAnalysis(classLabels = moduleColors,#基因所在的模块信息
                           identifiers = geneID$ENTREZID,
                           refCollection = GOcollection,
                           useBackground = "given",
                           threshold = 1e-4,
                           thresholdType = "Bonferroni",
                           getOverlapEntrez = TRUE,
                           getOverlapSymbols = TRUE,
                           ignoreLabels = "grey");
# 提取结果，并写入结果到文件
tab = GOenr$enrichmentTable
names(tab)
write.table(tab, file = "data/Step04-GOEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

# 提取主要结果，并写入文件
keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
# 小数位为2位
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# 给结果命名
colnames(screenTab) = c("module", "GOID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL

# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "data/Step04-GOEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

##### KEGG富集分析
# AnRichment没有直接提供KEGG数据的背景集
# 这里使用MSigDBCollection构建C2通路数据集
KEGGCollection = MSigDBCollection("data/msigdb_v7.1.xml", MSDBVersion = "7.1",
                                  organism = "rat",
                                  excludeCategories = c("h","C1","C3","C4","C5","C6","C7")) 
# KEGG分析
KEGGenr = enrichmentAnalysis(classLabels = moduleColors,
                           identifiers = geneID$ENTREZID,
                           refCollection = KEGGCollection,
                           useBackground = "given",
                           threshold = 1e-4,
                           thresholdType = "Bonferroni",
                           getOverlapEntrez = TRUE,
                           getOverlapSymbols = TRUE,
                           ignoreLabels = "grey")
# 提取KEGG结果，并写入文件
tab = KEGGenr$enrichmentTable
names(tab)
write.table(tab, file = "data/Step04-KEGGEnrichmentTable.xls", sep = "\t", quote = TRUE, row.names = FALSE)

# 提取主要结果并写入文件
keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
# 取两位有效数字
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# 对结果表格进行重命名
colnames(screenTab) = c("module", "ID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL
# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "data/Step04-KEGGEnrichmentTableScreen.xls", sep = "\t", quote = TRUE, row.names = FALSE)

### 输出cytoscape可视化
# 重新计算TOM，power值设置为前面选择好的
TOM = TOMsimilarityFromExpr(datExpr, power = 8)
# 输出全部网络模块
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = "data/Step04-CytoscapeInput-edges-all.txt",#基因间的共表达关系
                               nodeFile = "data/Step04-CytoscapeInput-nodes-all.txt",#
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = geneID$SYMBOL,
                               altNodeNames = geneID$ENTREZID,
                               nodeAttr = moduleColors)
# 输出感兴趣网络模块
modules = c("grey60","lightgreen","turquoise")
# 选择上面模块中包含的基因
inModule = is.finite(match(moduleColors, modules))
modGenes = geneID[inModule,]
# 选择指定模块的TOM矩阵
modTOM = TOM[inModule, inModule]

# 输出为Cytoscape软件可识别格式
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("data/Step04-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("data/Step04-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modGenes$SYMBOL,
                               altNodeNames = modGenes$ENTREZID,
                               nodeAttr = moduleColors[inModule])

##3、为了后续对指定模块中的基因进行下游分析，在此我们将每个模块的基因提取出来。
geneNames = colnames(datExpr)
cl=c()
moduleGenes=c()
module_colors=unique(moduleColors) 
for (color in module_colors){
     module=geneNames[which(moduleColors==color)]
     group1=rep(color,times=length(module))
     cl=c(cl,group1)
     moduleGenes=c(moduleGenes,module)
   }
my_module=data.frame("colour"=cl,"gene"=moduleGenes)
write.csv(my_module,"Modules_gene.csv")

