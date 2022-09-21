###1.准备工作####
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
#数据录入
expr <-read.delim("svminput.txt", row.names=NULL)
#创建任务
# 字符型转化为因子型
character2factor = function(x) {
  x %>%
    mutate(across(where(is.character),as.factor))
}
expr = character2factor(expr)
#创建任务
task = as_task_classif(expr,target = 'outcome', id = 'mdd', positive = 'treat')
print(task)
#创建学习器
learner <- lrn("classif.kknn", predict_type = "prob")
#定义内部循环重抽样
set.seed(123)
resampling <- rsmp("cv",folds = 10)
#定义模型评价指标
measure <- msr("classif.auc")
#搜索空间
search_space <- ps(k = p_int(lower = 1,upper = 10) )
#设定调优停止阈值
terminatior <- trm("none")
#设置调优算法
tuner <- tnr("random_search", batch_size = 20)
##打包调优需要的参数
at = AutoTuner$new(learner, resampling, measure, terminatior, tuner, search_space)#打包调优需要的参数
#外部重抽样
set.seed(123)
outer_resampling = rsmp("cv", folds = 5)
rr = resample(task, at, outer_resampling, store_models = TRUE)
#获取外部交叉验证结果（外部循环性能聚合）
rr$aggregate(measure)
extract_inner_tuning_results(rr)#获取超参数内部循环的表现（迭代次数假设为5次，那么就是在每一个外部循环的测试集上有5次超参数组合的测试，然后选择其中表现最优的超参数组合作为这个测试集代表性的超参数组合）
#这里其实可以通过示例数据可以有更好的理解

rr$score(measure)#获取超参数在外部循环的表现
rr$aggregate(msr("classif.auc"))#获得模型的准确率
rr$score(msr("classif.auc"))
rr1 <- rr$score(msr("classif.auc"))
write_excel_csv(rr1,file = 'rr_knn.csv')








