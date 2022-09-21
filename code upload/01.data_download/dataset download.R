rm(list=ls())
options(stringsAsFactors = F)##避免将character转换为因子
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
if(!require("GEOquery")) BiocManager::install("GEOquery")
if(!require("tidyverse")) BiocManager::install("tidyverse")
library(GEOquery)
library(dplyr)

gse <- getGEO("GSE98793", GSEMatrix = TRUE) 
show(gse)
class(gse)
str(gse)
a<-gse[[1]]
class(gse[[1]])##ExpressionSet

##提取第一个数据集的phenodata
dim(pData(gse[[1]]))
metdata<-pData(gse[[1]])
metdata[1:5,1:5]
colnames(metdata)##phenodata信息很多，但用得上的很少

##提取第一个表达矩阵
expma<-exprs(a)
dim(expma)
expma[1:5,1:5]
save(metdata,expma,file = "expma.Rdata")

load(file="expma.Rdata")

GPL="GPL570"##下载平台注释
gpl <- getGEO("GPL570",destdir = getwd()) 
GPL570 <- gpl@dataTable@table
save(GPL570, file = "GPL570.Rdata")
load(file = "GPL570.Rdata")
head(GPL570)
colnames(GPL570)

##取出注释信息
probe<-GPL570 %>% 
  select("ID","Gene Symbol","ENTREZ_GENE_ID")
head(probe)  
dim(probe)##22283个

library(tidyverse)
probe2<-apply(probe,1,function(x){
  paste(x[1],
        str_split(x[2]," /// ",simplify = T),##分割
        sep = "|")##连接符号
}  
) %>% 
  unlist()##得到的是个list
head(probe2)
length(probe2)##展开后得到24807个探针及对应关系
probe2<-as_tibble(probe2)

##注意这里的\\是用于转义 匹配分割 "|"，达到分割目的
probe2<-probe2 %>% separate(value,c("ID","GENE_SYMBOL"),sep = "\\|")
dim(probe2)##增加到24807行

## 找出重复ID，两个table的妙用
table(table(probe2$ID))##探针找出对应一个基因的有20878个，与grepl法得出的结果相同

## 下一步的目的即筛选出对应一个基因的探针
test2<-probe2 %>% count(ID) %>% filter(n==1) %>% ## count计数有点类似于table
  inner_join(probe2,by="ID") %>% ## 内连接只保留x,y中观测相同的变量
  select(-n)##remove "n" column
dim(test2)
head(test2)
probe2<-test2##将最终得到的结果赋值给probe2

write.csv(expma,file = 'GSE98793_dataframe.csv')
write.csv(probe2,file = 'GSE98793_id.csv')
write.csv(metdata,file = 'metadata.csv')
