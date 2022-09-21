#载入数据
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
##

if(T){
  load("expma.Rdata")
  load("Gpl570.Rdata")
}
expma[1:5,1:5]
expma1 <- expma
boxplot(expma)##看下表达情况

#metadata信息
metdata[1:5,1:5]


#探针注释信息
probe <- GPL570
head(probe)


#查看gene symbol是否有重复
table(duplicated(probe$`Gene Symbol`))


#整合注释信息到表达矩阵
ID<-rownames(expma)
expma<-apply(expma,1,function(x){log2(x+1)})
expma<-t(expma)
eset<-expma2[ID %in% probe$ID,] %>% cbind(probe)
eset[1:5,1:5]
expma2= 2^expma-1##反log2
colnames(eset)


#方法一：aggregate函数
test2<-aggregate(x=eset,by=list(eset$`Gene Symbol`),FUN=max,na.rm=T) 
test2##与去重结果相吻合

dim(test2)

colnames(test2)

write_csv(test2,file = '预处理数据.csv')

#方法二：dplyr
eset[1:5,1:5]

dim(eset)

colnames(eset)[8]<-"Gene"
colnames(eset)

test2<-eset[,c(1:6,8)] %>% 
  arrange(Gene) %>% 
  group_by(Gene) %>% 
  summarise_all(mean)
dim(test2)