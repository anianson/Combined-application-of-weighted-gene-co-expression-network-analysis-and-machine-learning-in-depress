library(tidyverse)
# 加载转录组数据
load(file = "All_data.Rdata")
colnames(gene_exp_res)
gene_exp_Ace <- gene_exp_res %>% 
  dplyr::select(starts_with("AceP") | starts_with("NC"))
sample_Ace <- sample_info %>% 
  dplyr::filter(group_res == "AceP" | group_res == "Control")
head(gene_exp_Ace)
head(sample_Ace)

# 样本添加batch信息
sample_Ace <- sample_Ace %>% 
  mutate(batch = paste0("batch",c(1,2,3,1,2,3)))
sample_Ace

# 第一张图：PCA图  fviz_pca_ind函数
expr_pca <- t(gene_exp_Ace)
pca_data <- cbind(expr_pca,sample_Ace)
expr_pca <- prcomp(expr_pca, scale = TRUE)
head(expr_pca)

library(factoextra)
# 按分组标注
fviz_pca_ind(expr_pca, 
             label="none", 
             habillage= pca_data$group_res,       # 按分组标注
             addEllipses=TRUE, 
             ellipse.level=0.95)
# 按批次标注
fviz_pca_ind(expr_pca, 
             label="none", 
             habillage= pca_data$batch,           # 按批次标注
             addEllipses=TRUE, 
             ellipse.level=0.95)

# 第二张图：聚类树图
# 按样本聚类
cluster_data <- t(gene_exp_Ace)
hc = hclust(dist(cluster_data))
plot(hc)

# 按batch聚类
cluster_data <- t(gene_exp_Ace)
rownames(cluster_data) = sample_Ace$batch
hc = hclust(dist(cluster_data))
plot(hc)

# 第三张图：聚类热图
# 列注释信息
annotation_col <- data.frame(sample_Ace$group_res,sample_Ace$batch)
rownames(annotation_col) <- colnames(gene_exp_Ace)

# 绘制热图
cg <- names(tail(sort(apply(gene_exp_Ace,1,sd)),1000))
library(pheatmap)
pheatmap(gene_exp_Ace[cg,],
         show_colnames = F, show_rownames = F,
         border_color = NA,
         annotation_col = annotation_col,
         # annotation_names_col = F,
         scale = "row")

# 去除批次效应  ComBat
# ??ComBat
library(sva)
mod <- model.matrix(~sample_Ace$group)
exprAce_batch <- ComBat(dat = gene_exp_Ace,
                        batch = sample_Ace$batch.ch1,
                        mod = mod)

# 第一张图：PCA图  fviz_pca_ind函数
expr_pca <- t(exprAce_batch)
pca_data <- cbind(expr_pca,sample_Ace)
expr_pca <- prcomp(expr_pca, scale = TRUE)
head(expr_pca)

library(factoextra)
# 按分组标注
fviz_pca_ind(expr_pca, 
             label="none", 
             habillage= pca_data$group,       # 按分组标注
             addEllipses=TRUE, 
             ellipse.level=0.95)
# 按批次标注
fviz_pca_ind(expr_pca, 
             label="none", 
             habillage= pca_data$batch.ch1,           # 按批次标注
             addEllipses=TRUE, 
             ellipse.level=0.95)

# 第二张图：聚类树图
# 按样本聚类
cluster_data <- t(exprAce_batch)
hc = hclust(dist(cluster_data))
plot(hc)

# 按batch聚类
cluster_data <- t(exprAce_batch)
rownames(cluster_data) = sample_Ace$batch.ch1
hc = hclust(dist(cluster_data))
plot(hc)

# 第三张图：聚类热图
# 列注释信息
annotation_col <- data.frame(sample_Ace$group,sample_Ace$batch.ch1)
rownames(annotation_col) <- colnames(exprAce_batch)

# 绘制热图
cg <- names(tail(sort(apply(exprAce_batch,1,sd)),1000))
library(pheatmap)
pheatmap(exprAce_batch[cg,],
         show_colnames = F, show_rownames = F,
         border_color = NA,
         annotation_col = annotation_col,
         # annotation_names_col = F,
         scale = "row")

# 小结：汇总一下去除批次效应前后的3张图，大家细品~！

write.csv(exprAce_batch, file='批次效应output.csv')
