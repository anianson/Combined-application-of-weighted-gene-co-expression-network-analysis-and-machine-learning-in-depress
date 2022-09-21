library(mlr3)
library(tidyverse)
library(kernlab)#提供内置数据集
library(mlr3measures)
library(mlr3learners)#提供学习器
library("mlr3tuning")#提供超参数调优
#数据录入
expr <-read.delim("svminput.txt", row.names=NULL)
str(expr)
character2factor = function(x) {
  x %>%
    mutate(across(where(is.character),as.factor))
}
expr = character2factor(expr)
expr$outcome <- as.factor(as.character(expr$outcome))
#创建任务
task <- as_task_classif(expr, target = "outcome")#创建任务
svm = lrn("classif.nnet",predict_type = "prob")#创建学习器
set.seed(123)
resampling <- rsmp("cv",folds = 10)
measure <- msr("classif.auc")
svm$param_set
svm$param_set$default
svm$param_set$levels
search_space = ps(
  maxit = p_int(lower = 1, upper = 500),
  size = p_int(lower = 1, upper = 10))

terminatior <- trm("evals", n_evals = 100)
tuner <- tnr("random_search")
at = AutoTuner$new(svm, resampling, measure, terminatior, tuner, search_space)
set.seed(123)
outer_resampling = rsmp("cv", folds = 5)
rr = resample(task, at, outer_resampling, store_models = TRUE)
rr$aggregate()
rr$aggregate(msr("classif.auc"))#获得模型的准确率
rr$score(msr("classif.auc"))
rr1 <- rr$score(msr("classif.auc"))
write_excel_csv(rr1,file = 'rr_nnet.csv')
