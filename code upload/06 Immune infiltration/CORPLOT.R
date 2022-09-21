library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggcorrplot)


data <- read.delim("~/Desktop/杨舒实验结果/A02-2022_07_09公共数据集生信/lassogene/环形热图inout.txt")
rownames(data) <- data[,1]
data <- data[,-1]

cordata <- round(cor(data), 3)
######
library(tidyverse)
library(reshape2)
df <- as.data.frame(cordata) %>% #将矩阵转换成数据框
  mutate(x=rownames(cordata)) %>%  #新建一列x，是11种属性变量
  melt(id='x') %>%                   #将宽数据转换成长数据，更适合ggplot2绘图
  rename('y'='variable','Corr'='value')  #将variable列名改成y

#先确定好绘图时x轴、y轴的因子顺序：
list <- rownames(cordata)
list <- factor(list,levels = list)
ggplot(df,aes(factor(x,levels = list),
                factor(y,levels = list), #定义x，y轴顺序，防止被默认改变
                fill=Corr))+  #根据相关性值填充颜色
  geom_tile()+  #色块函数
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red',
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL)+
  theme_bw(base_size = 15)
######
library(ggcorrplot)
ggcorrplot(cordata,lab=T)







                 