# https://bgreenwell.github.io/fastshap/articles/fastshap.html
setwd("C:/D/R project/Multi-Omics-Routine/Biostatistic/")

#设置种子数
set.seed(1011)

# 读取数据----------------------------------------------------------------------
load("myncRNA.Rdata")

sele <- c(sample(rownames(dat.raw)[dat.raw$group == 0], size = 2822, replace = FALSE),
          sample(rownames(dat.raw)[dat.raw$group == 1], size = 798, replace = FALSE))

dat.find <- dat.raw[sele,]
dat.vali <- dat.raw[!rownames(dat.raw) %in% sele, ]
dat.find <- dat.find %>% mutate(group = factor(group, labels = c("NC","CRC")))
dat.vali <- dat.vali %>% mutate(group = factor(group, labels = c("NC","CRC")))

# Step 3: 在测试集中，评估不同模型的预测能力，选择最佳的模型用于后续分析--------
# install.packages("shapviz")
# install.packages("shapley") # 目前的理解，只支持h2o包
# install.packages("fastshap")
# install.packages("iml")
# suppressMessages(require(shapley))
suppressMessages(require(shapviz))
suppressMessages(require(shapper))
suppressMessages(require(fastshap))
suppressMessages(require(iml))

model.lr <- glm(group ~ ., data = dat.find, family = binomial())
model.nb <- klaR::NaiveBayes(group ~ ., data = dat.find, fL = 0, usekernel = FALSE, adjust = 1)
group <- ifelse(dat.find$group == "NC", 0, 1)
dat4train <- dat.find[,1:2] %>% as.matrix
model.xgboost <- xgboost::xgboost(data = dat4train, label = group,
                                  max.depth = 3,
                                  eta = 0.1, 
                                  gamma = 1,
                                  nrounds = 1000, 
                                  min_child_weight = 0.5,
                                  colsample_bytree=1, subsmple = 1,
                                  nthread = 4,
                                  objective = "binary:logistic" # 二分类问题
)
# model.knn
model.rf <- randomForest::randomForest(group ~ ., data=dat.find, ntree=500, mtry=1, importance=FALSE, proximity=FALSE)
model.svm<- kernlab::ksvm(group ~ ., data = dat.find, type = "C-svc",
                          kernel="rbfdot", # Radial Basis kernel "Gaussian"
                          kpar=list(sigma=0.01), prob.model = TRUE, 
                          C = 10, na.action = na.omit, scaled = FALSE)
model.nn <- nnet::nnet(group ~ ., data = dat.find, size = 3, rang = 0.5,
                       decay = 0.01, maxit = 1000, linout = FALSE)
model.dt <- rpart::rpart(group ~ ., data = dat.find, control = rpart::rpart.control(cp = 0.0013), method = "class")

modelList <- list(
  model.lr, model.nb, model.xgboost, model.rf, model.svm, model.nn, model.dt
)
save(modelList, file = "Flow/modelList.Rdata")

# SHAP----
# 参考文献: 
# - 40567347
# SHAP分析应该用于测试集

# - kernelshap和fastshap结果基本一致
# - shapviz仅支持xgboost，lightgbm和h2o，xgboost计算结果和kernelshap，fastshap差别较大
# - shapper基于python运行，试运行后，程序卡死，无法输出结果
if(kernelshap){
  X_data <- dat.find[,1:2]
  Y_data <- dat.vali[,1:2]
  # LR----
  shap.lr <- kernelshap::kernelshap(model.lr, X = X_data) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.lr)
  # sv_dependence2D(shap.lr)
  sv_importance(shap.lr)
  sv_importance(shap.lr, kind = "beeswarm")
  # sv_interaction(shap.lr)
  sv_force(shap.lr)
  sv_waterfall(shap.lr) 
  # NB----
  shap.nb <- kernelshap::kernelshap(
    model.nb, X = X_data, 
    pred_fun = function(object, newdata){ predict(object, newdata)$posterior[,2] %>% as.numeric}) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.nb)
  # sv_dependence2D(shap.nb)
  sv_importance(shap.nb)
  sv_importance(shap.nb, kind = "beeswarm")
  # sv_interaction(shap.nb)
  sv_force(shap.nb)
  sv_waterfall(shap.nb) 
  # RF----
  shap.rf <- kernelshap::kernelshap(
    model.rf, X = X_data, 
    pred_fun = function(object, newdata){ predict(object, newdata, type = "prob")[,2] %>% as.numeric}) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.rf)
  # sv_dependence2D(shap.rf)
  sv_importance(shap.rf)
  sv_importance(shap.rf, kind = "beeswarm")
  # sv_interaction(shap.rf)
  sv_force(shap.rf)
  sv_waterfall(shap.rf) 
  # XGBoost----
  shap.xgb <- kernelshap::kernelshap(model.xgboost, X = as.matrix(X_data)) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.xgb)
  # sv_dependence2D(shap.xgb)
  sv_importance(shap.xgb)
  sv_importance(shap.xgb, kind = "beeswarm")
  # sv_interaction(shap.xgb)
  sv_force(shap.xgb)
  sv_waterfall(shap.xgb) 
  # SVM----
  shap.svm <- kernelshap::kernelshap(
    model.svm, X = X_data, 
    pred_fun = function(object, newdata){ predict(object, newdata, type = "probabilities")[,2] %>% as.numeric}) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.svm)
  # sv_dependence2D(shap.svm)
  sv_importance(shap.svm)
  sv_importance(shap.svm, kind = "beeswarm")
  # sv_interaction(shap.svm)
  sv_force(shap.svm)
  sv_waterfall(shap.svm) 
  # NN----
  shap.nn <- kernelshap::kernelshap(
    model.nn, X = X_data, 
    pred_fun = function(object, newdata){ predict(object, newdata)[,1] %>% as.numeric}) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.nn)
  # sv_dependence2D(shap.nn)
  sv_importance(shap.nn)
  sv_importance(shap.nn, kind = "beeswarm")
  # sv_interaction(shap.nn)
  sv_force(shap.nn)
  sv_waterfall(shap.nn) 
  # DT----
  shap.dt <- kernelshap::kernelshap(
    model.dt, X = X_data, 
    pred_fun = function(object, newdata){ predict(object, newdata)[,2] %>% as.numeric}) %>%
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.dt)
  # sv_dependence2D(shap.dt)
  sv_importance(shap.dt)
  sv_importance(shap.dt, kind = "beeswarm")
  # sv_interaction(shap.dt)
  sv_force(shap.dt)
  sv_waterfall(shap.dt) 
  # KNN---- 
}
if(shapper){
  X_data <- dat.find[,1:2]
  Y_data <- dat.vali[,1:2]
  # LR----
  shap.lr.py <- shapper::individual_variable_effect(
    model.lr, new_observation = X_data, data = X_data[1:100,]) %>% 
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.lr)
  # sv_dependence2D(shap.lr)
  sv_importance(shap.lr)
  sv_importance(shap.lr, kind = "beeswarm")
  # sv_interaction(shap.lr)
  sv_force(shap.lr)
  sv_waterfall(shap.lr) 
  
}
if(shapley){
  # only for h2o models
}
if(fastshap){
  # https://bgreenwell.github.io/fastshap/articles/fastshap.html
  X_data <- dat.find[,1:2]
  Y_data <- dat.vali[,1:2]
  # LR----
  shap.lr <- fastshap::explain(model.lr, X = X_data, newdata = X_data, 
                               adjust = TRUE, nsim = 100, shap_only=TRUE,
                               pred_wrapper = function(object, newdata){ predict(object, newdata, type = "response") }
  ) %>% 
    shapviz(object = ., X = X_data)
  # sv_dependence(shap.lr)
  # sv_dependence2D(shap.lr)
  sv_importance(shap.lr)
  sv_importance(shap.lr, kind = "beeswarm")
  # sv_interaction(shap.lr)
  sv_force(shap.lr)
  sv_waterfall(shap.lr) # Creates a waterfall plot of SHAP values of one observation. If multiple observations are selected, their SHAP values and predictions are averaged.
  
  # NB----
  shap.nb <- fastshap::explain(model.nb, X = X_data, newdata = Y_data, 
                               adjust = TRUE, nsim = 100, shap_only=TRUE,
                               pred_wrapper = function(object, newdata){ predict(object, newdata)$posterior[,2] %>% as.numeric}
  ) %>% 
    shapviz(object = ., X = Y_data)
  sv_importance(shap.nb)
  sv_importance(shap.nb, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.nb)
  sv_waterfall(shap.nb)
  # RF----
  shap.rf.fs <- fastshap::explain(model.rf, X = X_data, newdata = X_data,
                                  adjust = TRUE, nsim = 100, shap_only=TRUE,
                                  pred_wrapper = function(object, newdata){ predict(object, newdata, type = "prob")[,2] %>% as.numeric}
  ) %>% 
    shapviz(object = ., X = X_data)
  sv_importance(shap.rf.fs) 
  sv_importance(shap.rf.fs, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.rf.fs)
  sv_waterfall(shap.rf.fs)
  # XGBoost----
  shap.xgboost <- fastshap::explain(model.xgboost, X = as.matrix(X_data), newdata = as.matrix(Y_data),
                                    adjust = TRUE, nsim = 100, shap_only=TRUE,
                                    pred_wrapper = function(object, newdata){ predict(object, newdata)}
  ) %>% 
    shapviz(object = ., X = Y_data)
  sv_importance(shap.xgboost)
  sv_importance(shap.xgboost, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.xgboost)
  sv_waterfall(shap.xgboost)
  # SVM----
  shap.svm.fs <- fastshap::explain(model.svm, X = X_data, newdata = Y_data,
                                   adjust = TRUE, nsim = 100, shap_only=TRUE,
                                   pred_wrapper = function(object, newdata){ predict(object, newdata, type = "probabilities")[,2] %>% as.numeric}
  ) %>% 
    shapviz(object = ., X = Y_data)
  sv_importance(shap.svm.fs)
  sv_importance(shap.svm.fs, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.svm.fs)
  sv_waterfall(shap.svm.fs)
  # NN----
  shap.nn <- fastshap::explain(model.nn, X = X_data, newdata = Y_data,
                               adjust = TRUE, nsim = 100, shap_only=TRUE,
                               pred_wrapper = function(object, newdata){ predict(object, newdata)[,1] %>% as.numeric}
  ) %>% 
    shapviz(object = ., X = X_data)
  sv_importance(shap.nn)
  sv_importance(shap.nn, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.nn)
  sv_waterfall(shap.nn)
  # DT----
  shap.dt <- fastshap::explain(model.dt, X = X_data, newdata = Y_data,
                               adjust = TRUE, nsim = 100, shap_only=TRUE,
                               pred_wrapper = function(object, newdata){ predict(object, newdata)[,2] %>% as.numeric}
  ) %>% 
    shapviz(object = ., X = Y_data)
  sv_importance(shap.dt)
  sv_importance(shap.dt, kind = "beeswarm") # kind = c("bar", "beeswarm", "both", "no")
  sv_force(shap.dt)
  sv_waterfall(shap.dt)
  # KNN----
}
if(shapviz){
  # only for XGBoost, LightGBM or H2O
  myshap <- shapviz::shapviz(model.xgboost, X_pred = dat.vali[,1:2] %>% as.matrix)
  
  sv_dependence(myshap, c("hsa_miR_744_5p","hsa_miR_1247_3p"))
  
  sv_importance(myshap)
  # sv_importance(shap.xgboost)
  
  sv_importance(myshap, kind = "beeswarm")
  sv_importance(shap.xgboost, kind = "beeswarm")
  
  sv_force(myshap)
  sv_force(shap.xgboost)
  
  sv_waterfall(myshap)
  sv_waterfall(shap.xgboost)
}













