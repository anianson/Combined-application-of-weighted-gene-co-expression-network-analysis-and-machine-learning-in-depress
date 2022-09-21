library('ComplexHeatmap')

library('circlize')
##读入数据
环形热图inout <- read.delim("~/Desktop/杨舒实验结果/A02-2022_07_09公共数据集生信/lassogene/环形热图inout.txt")
data <- 环形热图inout
head(data)
row.names(data)<-data[,1]#修改行名

data<-data[,-1]##删除第一列，使之变为数字矩阵，绘图的数据要求为矩阵（也就是单一类型的数据矩阵，这里全为数字）
##转化为矩阵并对其进行归一化

madt<-as.matrix(data)

madt2<-t(scale(t(madt)))

Heatmap(madt2)     ##这就是默认参数，普通的热图

range(madt2)##查看值的分布

mycol=colorRamp2(c(-2, 0.15, 2.3),c("blue", "white", "red")) ##定义颜色范围

circos.heatmap(madt2,col=mycol)##默认参数

circos.clear()  ##绘制完成后需要使用此函数完全清除布局

circos.par(gap.after=c(50))    ##调整圆环首尾间的距离，数值越大，距离越宽

circos.heatmap(madt2,col=mycol, dend.side="inside", cluster=TRUE)


library(dendextend)

library(dendsort)

circos.par(gap.after=c(50))

circos.heatmap(madt2,col=mycol,dend.side="inside",
               track.height=0.38,cluster=TRUE,
               dend.track.height=0.18,
               dend.callback=function(dend,m,si){color_branches(dend,k=15,col=1:15)})

##添加图例标签等

lg=Legend(title="Exp",col_fun=mycol,direction= c("vertical"))

grid.draw(lg)


column_od = hclust(dist(t(madt2)))$order #对列聚类
circos.track(track.index = get.current.track.index(), ##将列名添加在第二个轨道（就是热图所在的环形轨道）
             
             panel.fun = function(x, y) {
               
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 
                 cn = colnames(madt2)##取得列名
                 
                 n = length(cn)
                 
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.3, "mm"), ##x轴坐标
                             
                             4.1+(1:n)*0.7, ##y轴坐标
                             
                             cn, ##输入要展示的列名
                             
                             cex = 0.3, ##列名的大小
                             
                             adj = c(0, 1),
                             
                             facing = "inside")
                 
               }
               
             }, bg.border = NA)

circos.clear()


