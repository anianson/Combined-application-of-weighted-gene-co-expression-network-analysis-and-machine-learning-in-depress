library(tidyverse)
library(ggplot2)

data <- read.table("/A02-2022_07_09公共数据集生信/4.差异表达改/all.txt", header = T)

data$col <- "no significant"
data$col[data$adj.P.Val < 0.05 & data$logFC > 0] <- "Up"
data$col[data$adj.P.Val < 0.05 & data$logFC < 0] <- "Down"
data$col <- factor(data$col, levels = c("Down", "no significant","Up"))

data$size <- 1

ggplot() +
  geom_point(data = data, aes(logFC, -log10(adj.P.Val), colour = col, fill = col),
             size = data$size) +
  scale_colour_manual(values = c("#4DBBD5", "grey", "#E64B35"))  +
  geom_vline(xintercept = 0, color="grey40", linetype=2) +
  geom_hline(yintercept = -log10(0.05), color="grey40", linetype=2)
