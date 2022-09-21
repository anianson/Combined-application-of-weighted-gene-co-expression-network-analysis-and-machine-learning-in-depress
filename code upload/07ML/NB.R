#贝叶斯模型
####1.准备工作####
library(mlr3)#机器学习R包
library(mlr3learners)
library(tidyverse)#数据处理R包
library(mlr3verse)
#数据录入
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
learner <- lrn("classif.naive_bayes", predict_type = "prob")
set.seed(123)
resampling <- rsmp("cv", folds=10)#默认为重复一次的10折交叉验证
rr = resample(task,learner, resampling, store_models = TRUE)
rr$aggregate(msr("classif.auc"))#获得模型的准确率
rr$score(msr("classif.auc"))
rr1 <- rr$score(msr("classif.auc"))
write_excel_csv(rr1,file = 'rr_NB.csv')







