library(reshape2)
library(ggpubr)
File='多组箱线图input.txt'
rt=read.table(File, header=T, sep="\t", check.names=F, row.names=1)
rt$sample_type=factor(rt$sample_type, levels=c("con","treat"))
colnames(rt)

data=melt(rt, id.vars=c("sample_type"))
colnames(data)=c("sample_type", "gene_ID", "expression_level")
p=ggboxplot(data, x="gene_ID", y="expression_level", fill = "sample_type",
            xlab="",
            ylab="expression_level",
            legend.title="group",
            palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
show(p)
p1=p+stat_compare_means(aes(group=sample_type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
show(p1)
#
pdf(file="boxplot1.pdf", width=6, height=5)
print(p1)
dev.off()
