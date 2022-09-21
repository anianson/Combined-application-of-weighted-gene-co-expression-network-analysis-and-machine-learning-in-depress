library(mlr3)
library(mlr3verse)
library(mlr3viz)
library(mlr3learners)
library(mlr3tuning)
library(tidyverse)
#导数据
expr <-read.delim("svminput.txt", row.names=NULL)
# 字符型转化为因子型
character2factor = function(x) {
  x %>%
    mutate(across(where(is.character),as.factor))
}
expr = character2factor(expr)
task <- as_task_classif(expr, target = "outcome")#创建任务
learner <- lrn("classif.ranger", predict_type = "prob")#创建学习器
set.seed(123)
resampling <- rsmp("cv",folds = 10)#定义内部循环重抽样

measure <- msr("classif.auc")#定义模型评价指标
search_space = ps(
  num.trees = p_int(lower = 300, upper = 300),
  mtry = p_int(lower = 6, upper = 12),
  min.node.size = p_int(lower = 1,upper = 5),
  max.depth = p_int(lower = 5,upper = 20)
)#定义超参数空间

evals1000 = trm("evals", n_evals = 100)#设定调优停止阈值

tuner = tnr("random_search")#设置调优算法

at = AutoTuner$new(learner, resampling, measure, evals1000 , tuner, search_space)#打包调优需要的参数
#这里打包的是基于内部循环的，因为我们的超参数调优就是内部循环的工作
set.seed(123)
outer_resampling = rsmp("cv", folds = 5)

rr = resample(task, at, outer_resampling, store_models = TRUE)

rr$aggregate()#获取外部交叉验证结果（外部循环性能聚合）
rr$aggregate(msr("classif.auc"))#获得模型的准确率
rr$score(msr("classif.auc"))
rr1 <- rr$score(msr("classif.auc"))
write_excel_csv(rr1,file = 'rr_random.csv')









