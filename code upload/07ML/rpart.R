library(tidyverse)
library(mlr3)
library(mlr3verse)
library(mlr3learners)
library(rpart.plot)
#导数据
expr <-read.delim("svminput.txt", row.names=NULL)
# 字符型转化为因子型
character2factor = function(x) {
  x %>%
    mutate(across(where(is.character),as.factor))
}
expr = character2factor(expr)
task <- as_task_classif(expr, target = "outcome")#创建任务
learner <- lrn("classif.rpart", predict_type = "prob")#创建学习器
set.seed(123)
resampling <- rsmp("cv",folds = 10)#定义内部循环重抽样

measure <- msr("classif.specificity")#定义模型评价指标
measures = msrs(c("classif.auc","classif.acc","classif.ce","classif.sensitivity","classif.specificity")) 
search_space = ps(
  cp = p_dbl(lower = 0.001, upper = 0.1),
  minsplit = p_int(lower = 1, upper = 10),
  maxdepth = p_int(lower = 3,upper = 10),
  minbucket = p_int(lower = 3,upper = 10)
)#定义超参数空间

evals1000 = trm("evals", n_evals = 100)#设定调优停止阈值

tuner = tnr("random_search")#设置调优算法

at = AutoTuner$new(learner, resampling, measure, evals1000 , tuner, search_space)#打包调优需要的参数
#这里打包的是基于内部循环的，因为我们的超参数调优就是内部循环的工作
set.seed(123)
outer_resampling = rsmp("cv", folds = 5)

rr = resample(task, at, outer_resampling, store_models = TRUE)
#嵌套重采样并不局限于超参数调优。您可以将AutoTuner替换为AutoFSelector，并估计拟合在优化的特征子集上的模型的性能
#时间大概在20min左右

rr$aggregate()#获取外部交叉验证结果（外部循环性能聚合）
rr$aggregate(msr("classif.specificity"))#获得模型的准确率
rr$score(msr("classif.specificity"))
rr1 <- rr$score(msr("classif.specificity"))
write_excel_csv(rr1,file = 'rr_rpart-classif.specificity.csv')









