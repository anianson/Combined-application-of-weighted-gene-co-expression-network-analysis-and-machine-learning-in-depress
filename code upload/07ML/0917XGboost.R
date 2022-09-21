library(tidyverse)
library(mlr3)
#https://mlr3book.mlr-org.com/basics.html#tasks
library(tidymodels)
library(ggplot2)
library("mlr3verse")
library(kknn)
library(paradox)#定义超参数的搜索空间
library(mlr3tuning)#用于对超参数调优
library(ggplot2)
expr <-read.delim("svminput.txt", row.names=NULL)#导入数据集

#XGboost预测变量需要是数值型变量
expr$outcome <- as.factor(as.character(expr$outcome))

task <- as_task_classif(expr, target = "outcome")#创建任务
lrn <- lrn("classif.xgboost",predict_type = "prob")

lrn$param_set

search_space = ps(
  eta = p_int(lower = 0, upper = 1),
  gamma = p_int(lower = 0, upper = 5),
  max_depth = p_int(lower = 1,upper = 5),
  nrounds = p_int(lower = 800,upper = 800),
  subsample = p_dbl(lower = 0.5,upper = 1),
  eval_metric = p_fct(c('auc', 'logloss'))
)
set.seed(123)
resampling <- rsmp("cv",folds = 10)#定义内部循环重抽样
measure <- msr("classif.auc")#定义模型评价指标
evals1000 = trm("evals", n_evals = 100)#设定调优停止阈值
tuner = tnr("random_search")#设置调优算法
at = AutoTuner$new(lrn, resampling, measure, evals1000 , tuner, search_space)#打包调优需要的参数
set.seed(123)
outer_resampling = rsmp("cv", folds = 5)
ptm <- proc.time()
rr = resample(task, at, outer_resampling, store_models = TRUE)
proc.time() - ptm
rr$aggregate(msr("classif.auc"))#获得模型的准确率
rr$score(msr("classif.auc"))
rr1 <- rr$score(msr("classif.auc"))
write_excel_csv(rr1,file = 'rr_XGboost.csv')





