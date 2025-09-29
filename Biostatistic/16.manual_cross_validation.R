# Step 1: 按照7:3的比例，将数据拆分为内部训练和内部测试集-----------------------
# install.packages("sampling") # sample
# install.packages("xgboost") # xgboost
# install.packages("gbm") # gbm
# install.packages("e1071") # NaiveBayes, SVM
# install.packages("kernlab") # SVM
# install.packages("klaR") # NaiveBayes
# install.packages("randomForest") # RandomForest
# install.packages("ranger") # RandomForest
# install.packages("nnet") # NNET
# install.packages("neuralnet") # NNET
# install.packages("tensorflow") # NNET
# install.packages("keras") # NNET
# install.packages("class") # KNN
# install.packages("FNN") # KNN
# install.packages("kknn") # KNN
# install.packages("rpart") # DT
# install.packages("caret") # cross validation & tune
# install.packages("tidymodels") # cross validation & tune

suppressMessages(require(dplyr))
suppressMessages(require(ggplot2))
suppressMessages(require(stringr))
suppressMessages(require(caret))

in_training <- createDataPartition(dat.raw$group, p = .7, list = FALSE)
training_data <- dat.raw[in_training, ]
testing_data  <- dat.raw[-in_training, ]

save(dat.raw, in_training, training_data, testing_data, file = "Flow/dat4model.Rdata")
# Step 2: 在内部训练集中，使用嵌套交叉验证（蒙特卡洛或10折），评估不同模型的预测能力----
# https://topepo.github.io/caret/available-models.html

set.seed(1011)

# 10 fold cv sampling-----------------------------------------------------------
# 尽管蒙特卡洛交叉验证是一个很好的选择，但是运算量太大了，所以，选择10 fold cv
# 0.9 * nrow(training_data) # [1] 4562.1
# 0.1 * nrow(training_data) # [1] 506.9

folds_list <- createFolds(
  y = training_data$group,
  k = 10,
  list = TRUE,
  returnTrain = FALSE
)

# lapply(folds_list, length)
save(folds_list, file = "Flow/foldlist.Rdata")

load("Flow/dat4model.Rdata")
load("Flow/foldlist.Rdata")
# 1. LR-------------------------------------------------------------------------
res.lr <- data.frame(
  Fold = paste0("Fold", c(1:10)),
  stringsAsFactors = FALSE
) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA)

res.lr[,-1] <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  fit.lr <- glm(group ~ ., data = in.train, family = binomial())
  modelroc.lr <- roc(in.test$group, predict(fit.lr, newdata = in.test[,1:2], type='response'))
  
  index <- which.max(modelroc.lr$sensitivities + modelroc.lr$specificities - 1)
  
  res <- c(modelroc.lr$auc, 
           modelroc.lr$sensitivities[index], 
           modelroc.lr$specificities[index])
  names(res) <- c("auc","sens","spec")
  return(res)
}) %>% as.data.frame %>% t

# control = trainControl(method = "cv", 
#                        number = 10, 
#                        preProcOptions = NULL, # list(method = c("center", "scale")),
#                        classProbs = TRUE, 
#                        summaryFunction = twoClassSummary,
#                        seed = c(rep(1011, 10), 1011),
#                        returnResamp = "all",
#                        allowParallel = FALSE)
# lr.train = train(group ~ ., data = training_data %>% 
#                    mutate(group = factor(group, labels = c("NC", "CRC"))),
#                  method = "glm", metric = "ROC",
#                  family = binomial(),
#                  trControl = control)
# lr.train$resample
# 2. {scale-insens} NB-------------------------------------------------------------------------
res.nb <- data.frame(
  Fold = paste0("Fold", c(1:10)),
  stringsAsFactors = FALSE
) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA)

res.nb[,-1] <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  fit.nb <- naiveBayes(group ~ ., data = in.train, laplace=0)
  modelroc.nb <- roc(in.test$group, predict(fit.nb, newdata = in.test[,1:2]) %>% as.numeric)
  
  index <- which.max(modelroc.nb$sensitivities + modelroc.nb$specificities - 1)
  
  res <- c(modelroc.nb$auc, 
           modelroc.nb$sensitivities[index], 
           modelroc.nb$specificities[index])
  names(res) <- c("auc","sens","spec")
  return(res)
}) %>% as.data.frame %>% t
# 3. {scale-sens} KNN [Tune]------------------------------------------------------------------------
res.knn <- data.frame(
  Fold = rep(paste0("Fold", c(1:10)), each = 10),
  para.k = rep(seq(5, 50, 5), 10),
  stringsAsFactors = FALSE
) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA)

res.knn[,3:5] <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  klist <- seq(5, 50, 5) %>% as.list
  names(klist) <- seq(5, 50, 5)
  lapply(klist, function(k){
    knn_pred <- knn(train = in.train[,1:2], test = in.test[,1:2],
                    cl = in.train$group, k = k)
    modelroc.knn <- roc(in.test$group, knn_pred %>% as.numeric)
    
    index <- which.max(modelroc.knn$sensitivities + modelroc.knn$specificities - 1)
    
    res <- c(modelroc.knn$auc, 
             modelroc.knn$sensitivities[index], 
             modelroc.knn$specificities[index])
    names(res) <- c("auc","sens","spec")
    return(res)
  })
}) %>% as.data.frame %>% t
# 4. {tree} RF [Tune]-------------------------------------------------------------------------
res.rf <- data.frame(
  Fold = rep(paste0("Fold", c(1:10)), each = 2),
  para.mtry = rep(c(1,2), 10),
  stringsAsFactors = FALSE
) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA)

res.rf[,3:5] <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  mtrylist <- list(1,2)
  names(mtrylist) <- mtrylist
  lapply(mtrylist, function(mtry){
    fit.rf <- randomForest(group ~ ., data=in.train, ntree=1000, mtry=mtry,
                           importance=FALSE, proximity=FALSE)
    modelroc.rf <- roc(in.test$group, predict(fit.rf, in.test) %>% as.numeric)
    
    index <- which.max(modelroc.rf$sensitivities + modelroc.rf$specificities - 1)
    
    res <- c(modelroc.rf$auc, 
             modelroc.rf$sensitivities[index], 
             modelroc.rf$specificities[index])
    names(res) <- c("auc","sens","spec")
    return(res)
  })
}) %>% as.data.frame %>% t
# 5. {scale-sens} SVM[Tune]------------------------------------------------------------------------
res.svm = expand.grid(Fold = paste0("Fold", c(1:10)),
                      para.cost = c(0.01,0.1,1,10,100),
                      para.gamma= c(0.01,0.1,1,10,100)) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA) %>% 
  arrange(Fold)

tmp.svm <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  costlist <- list(0.01,0.1,1,10,100)
  names(costlist) <- paste0("cost:",costlist)
  lapply(costlist, function(cost){
    
    gammalist <- list(0.01,0.1,1,10,100)
    names(gammalist) <- paste0("gamma:",gammalist)
    lapply(gammalist,function(gamma){
      fit.svm <- e1071::svm(group ~ ., data=in.train, cost=cost, gamma = gamma, kernel = "radial")
      modelroc.svm <- roc(in.test$group, predict(fit.svm, in.test[,1:2]) %>% as.numeric)
      
      index <- which.max(modelroc.svm$sensitivities + modelroc.svm$specificities - 1)
      
      res <- c(modelroc.svm$auc, 
               modelroc.svm$sensitivities[index], 
               modelroc.svm$specificities[index])
      names(res) <- c("auc","sens","spec")
      return(res)
    })
  })
}) %>% as.data.frame %>% t
res.svm[,4:6] <- tmp.svm

res.svm %>% 
  group_by(para.cost, para.gamma) %>% 
  mutate(mean = mean(Auc)) %>% 
  ggplot(aes(x = para.cost %>% factor, y = para.gamma %>% factor, fill = mean)) +
  geom_tile() +
  scale_fill_viridis_c()
res.svm %>% 
  group_by(para.cost, para.gamma) %>% 
  mutate(mean = mean(Auc)) %>% 
  arrange(desc(mean))
# 6. {tree} XGBoost[Tune]--------------------------------------------------------------------
res.xgb = expand.grid(Fold = paste0("Fold", c(1:10)),
                      max_depth = c(2:10),
                      eta = c(0.05, 0.1, 0.2, 0.3),
                      gamma = c(0.01, 0.05, 0.1, 0.2, 0.5, 1.0),
                      colsample_bytree = c(0.5, 1),
                      min_child_weight = c(0.5, 1, 1.5)) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA) %>% 
  arrange(min_child_weight) %>% 
  arrange(colsample_bytree) %>% 
  arrange(gamma) %>% 
  arrange(eta) %>% 
  arrange(max_depth) %>% 
  arrange(Fold)

tmp.xgb <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  max_depthlist <- as.list(c(2:10))
  names(max_depthlist) <- paste0("max_depth:", max_depthlist)
  lapply(max_depthlist, function(max_depth){
    
    etalist <- list(0.05, 0.1, 0.2, 0.3)
    names(etalist) <- paste0("eta:",etalist)
    lapply(etalist, function(eta){
      
      gammalist <- list(0.01, 0.05, 0.1, 0.2, 0.5, 1.0)
      names(gammalist) <- paste0("gamma:",gammalist)
      lapply(gammalist,function(gamma){
        
        colsample_bytreelist <- list(0.5, 1)
        names(colsample_bytreelist) <- paste0("colsample_bytree:",colsample_bytreelist)
        lapply(colsample_bytreelist,function(colsample_bytree){
          
          min_child_weightlist <- list(0.5, 1, 1.5)
          names(min_child_weightlist) <- paste0("min_child_weight:",min_child_weightlist)
          lapply(min_child_weightlist,function(min_child_weight){
            suppressMessages(fit.xgb <- xgboost(data = in.train[,1:2] %>% as.matrix, 
                                                label = in.train$group %>% as.character %>% as.numeric,
                                                subsample = 1, nrounds = 50, 
                                                colsample_bytree=colsample_bytree,
                                                max.depth = max_depth, # 通常 3 ~ 6
                                                eta = eta, # 通常 0.01 ~ 0.3
                                                gamma = gamma, min_child_weight = min_child_weight,
                                                nthread = 1,
                                                objective = "binary:logistic" # 二分类问题
            ))
            modelroc.xgb <- roc(in.test$group, predict(fit.xgb, in.test[,1:2] %>% as.matrix))
            
            index <- which.max(modelroc.xgb$sensitivities + modelroc.xgb$specificities - 1)
            
            res <- c(modelroc.xgb$auc, 
                     modelroc.xgb$sensitivities[index], 
                     modelroc.xgb$specificities[index])
            names(res) <- c("auc","sens","spec")
            return(res)
          })
        })
      })
    })
  })
}) %>% as.data.frame %>% t
res.xgb[,7:9] <- tmp.xgb

res.xgb %>% 
  group_by(max_depth, eta, gamma, colsample_bytree, min_child_weight) %>% 
  mutate(mean = mean(Auc)) %>% 
  # distinct(max_depth, eta, gamma, colsample_bytree, min_child_weight, .keep_all = TRUE) %>%
  arrange(desc(mean))
# 7. {tree} DT[Tune]-------------------------------------------------------------------------
res.dt <- expand.grid(
  Fold = paste0("Fold", c(1:10)),
  para.cp = seq(0.0001, .05, by = 0.0001),
  stringsAsFactors = FALSE
) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA)

res.dt[,3:5] <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  cplist <- as.list(seq(0.0001, .05, by = 0.0001))
  names(cplist) <- cplist
  lapply(cplist, function(cp){
    fit.dt <- rpart::rpart(group~., data=in.train, method = "class", 
                           control = rpart::rpart.control(cp = cp))
    modelroc.dt <- roc(in.test$group, predict(fit.dt, in.test[,1:2], type="class") %>% as.numeric)
    ci(modelroc.dt)
    
    index <- which.max(modelroc.dt$sensitivities + modelroc.dt$specificities - 1)
    
    res <- c(modelroc.dt$auc, 
             modelroc.dt$sensitivities[index], 
             modelroc.dt$specificities[index])
    names(res) <- c("auc","sens","spec")
    return(res)
  })
}) %>% as.data.frame %>% t
# 8. {scale-sens} NN-----------------------------------------------------------------------
res.nn = expand.grid(
  Fold = paste0('Fold', c(1:10)),
  size = c(1:20),
  decay = c(0, 1e-5, 1e-4, 1e-3, 1e-2)) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA) %>% 
  arrange(decay) %>% 
  arrange(size) %>% 
  arrange(Fold)

tmp.nn <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  sizelist <- as.list(c(1:20))
  names(sizelist) <- paste0("size:",sizelist)
  lapply(sizelist, function(size){
    
    decaylist <- as.list(c(0, 1e-5, 1e-4, 1e-3, 1e-2))
    names(decaylist) <- paste0("decay:",decaylist)
    lapply(decaylist, function(decay){
      fit.nn <- nnet::nnet(group ~ ., data = in.train, size = size, 
                           decay = decay, maxit = 1000, linout = FALSE)
      modelroc.nn <- roc(in.test$group, predict(fit.nn, in.test[,1:2])[,1] %>% as.numeric)
      
      index <- which.max(modelroc.nn$sensitivities + modelroc.nn$specificities - 1)
      
      res <- c(modelroc.nn$auc, 
               modelroc.nn$sensitivities[index], 
               modelroc.nn$specificities[index])
      names(res) <- c("auc","sens","spec")
      return(res)
    })
  })
}) %>% as.data.frame %>% t
res.nn[,4:6] <- tmp.nn
# * 9. {scale-sens} NNET-----------------------------------------------------------------------
# mulitiple layers
# fit.neun <- neuralnet::neuralnet(
#   group == 1 ~ ., dat.find, rep = 1, hidden = c(3,3),
#   linear.output = FALSE)
# modelroc.neun <- roc(dat.vali$group, predict(fit.neun, dat.vali[,1:2])[,1] %>% as.numeric)
# ci(modelroc.neun)
# plot(fit.neun)
# NeuralNetTools::plotnet(fit.neun)

res.neun = expand.grid(Fold = paste0("Fold", c(1:10)),
                       hidden = c(1:20)) %>% 
  mutate(Auc = NA, Sens = NA, Spec = NA) %>% 
  arrange(Fold)

tmp.neun <- lapply(folds_list, function(x){
  in.train <- training_data[-x,]
  in.test  <- training_data[x,]
  
  preProc <- preProcess(in.train[,1:2], method = c("center", "scale"))
  in.train[, 1:2]<- predict(preProc, in.train[, 1:2])
  in.test[, 1:2] <- predict(preProc, in.test [, 1:2])
  
  lapply(as.list(c(1:20)), function(hidden){
    fit.neun <- neuralnet::neuralnet(
      group == 1 ~ ., in.train, rep = 1, hidden = hidden,
      linear.output = FALSE)
    modelroc.neun <- roc(in.test$group, predict(fit.neun, in.test[,1:2])[,1] %>% as.numeric)
    
    index <- which.max(modelroc.neun$sensitivities + modelroc.neun$specificities - 1)
    
    res <- c(modelroc.neun$auc, 
             modelroc.neun$sensitivities[index], 
             modelroc.neun$specificities[index])
    names(res) <- c("auc","sens","spec")
    return(res)
  })
})
tmp.neun <- tmp.neun %>% as.data.frame %>% t
res.neun[,3:5] <- tmp.neun
# RES-------------------------------------------------------------------------
res <- list(LR = res.lr, 
            NB = res.nb,
            KNN = res.knn,
            RF = res.rf,
            SVM = res.svm,
            XGBoost = res.xgb,
            NN = res.nn,
            DT = res.dt
            # , NN1 = res.neun
)
save(res, file = "Flow/ML(scale).Rdata")

# Best model--------------------------------------------------------------------
# 1.LR----
# 2.NB----
# 3.KNN----
res.knn %>% 
  group_by(para.k) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(para.k)) %>% 
  arrange(desc(mean))
# # A tibble: 10 × 6
# Fold  para.k   Auc  Sens  Spec  mean
# <chr>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1     50 0.858 0.741 0.975 0.861
# 2 Fold1     40 0.858 0.741 0.975 0.860
# 3 Fold1     45 0.853 0.732 0.975 0.860
# 4 Fold1     35 0.864 0.75  0.977 0.859
# 5 Fold1     20 0.857 0.741 0.972 0.859
# 6 Fold1     15 0.843 0.714 0.972 0.859
# 7 Fold1     30 0.857 0.741 0.972 0.858
# 8 Fold1     25 0.849 0.723 0.975 0.858
# 9 Fold1     10 0.846 0.723 0.970 0.858
# 10 Fold1    5 0.843 0.732 0.954 0.857
# 4.RF----
res.rf %>% 
  group_by(para.mtry) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(para.mtry)) %>% 
  arrange(desc(mean))
# # A tibble: 2 × 6
# Fold  para.mtry   Auc  Sens  Spec  mean
# <chr>     <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1         1 0.871 0.777 0.965 0.860
# 2 Fold1         2 0.863 0.759 0.967 0.857
# 5.SVM----
res.svm %>% 
  group_by(para.cost, para.gamma) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(paste0(para.cost, para.gamma))) %>% 
  arrange(desc(mean))
# # A tibble: 25 × 7
# Fold  para.cost para.gamma   Auc  Sens  Spec  mean
# <fct>     <dbl>      <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1      1         100   0.867 0.759 0.975 0.860
# 2 Fold1     10          10   0.842 0.714 0.970 0.860
# 3 Fold1     10           1   0.852 0.732 0.972 0.859
# 4 Fold1     10         100   0.836 0.705 0.967 0.857
# 5 Fold1      1          10   0.862 0.75  0.975 0.857
# 6 Fold1      1           1   0.864 0.75  0.977 0.856
# 7 Fold1      1           0.1 0.860 0.741 0.980 0.854
# 8 Fold1      0.1       100   0.850 0.723 0.977 0.850
# 9 Fold1      0.01        1   0.849 0.723 0.975 0.850
# 10 Fold1      0.01      100   0.857 0.732 0.982 0.850
# # ℹ 15 more rows
# # ℹ Use `print(n = ...)` to see more rows
# 6.XGBoost----
res.xgb %>% 
  group_by(max_depth, eta, gamma, colsample_bytree, min_child_weight) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(paste0(max_depth, eta, gamma, colsample_bytree, min_child_weight))) %>% 
  arrange(desc(mean))
# # A tibble: 1,296 × 10
# Fold  max_depth   eta gamma colsample_bytree min_child_weight   Auc  Sens  Spec  mean
# <fct>     <int> <dbl> <dbl>            <dbl>            <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1         2   0.2  0.01              0.5              0.5 0.953 0.893 0.916 0.949
# 2 Fold1         2   0.2  1                 0.5              1   0.954 0.902 0.916 0.949
# 3 Fold1         3   0.1  0.01              1                1.5 0.954 0.884 0.919 0.949
# 4 Fold1         3   0.1  0.01              1                0.5 0.953 0.884 0.914 0.949
# 5 Fold1         3   0.1  0.01              1                1   0.954 0.884 0.914 0.949
# 6 Fold1         3   0.1  0.05              1                0.5 0.953 0.884 0.914 0.949
# 7 Fold1         3   0.1  0.5               1                0.5 0.953 0.884 0.911 0.949
# 8 Fold1         3   0.1  0.5               1                1.5 0.954 0.884 0.919 0.949
# 9 Fold1         3   0.1  1                 1                1.5 0.953 0.884 0.916 0.949
# 10 Fold1         3   0.1  0.05              1                1   0.954 0.884 0.914 0.949
# # ℹ 1,286 more rows
# # ℹ Use `print(n = ...)` to see more rows
# 7.DT----
res.dt %>% 
  group_by(para.cp) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(para.cp)) %>% 
  arrange(desc(mean))
# # A tibble: 500 × 6
# Fold  para.cp   Auc  Sens  Spec  mean
# <chr>   <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1  0.0405 0.890 0.811 0.970 0.890
# 2 Fold1  0.0406 0.890 0.811 0.970 0.890
# 3 Fold1  0.0407 0.890 0.811 0.970 0.890
# 4 Fold1  0.0408 0.890 0.811 0.970 0.885
# 5 Fold1  0.0403 0.884 0.802 0.967 0.884
# 6 Fold1  0.0351 0.884 0.804 0.965 0.883
# 7 Fold1  0.0404 0.881 0.793 0.970 0.882
# 8 Fold1  0.0402 0.884 0.802 0.967 0.882
# 9 Fold1  0.0409 0.881 0.793 0.970 0.881
# 10 Fold1  0.041  0.881 0.793 0.970 0.881
# # ℹ 490 more rows
# # ℹ Use `print(n = ...)` to see more rows
# 8.NN----
res.nn %>% 
  group_by(size, decay) %>% 
  mutate(mean = mean(Auc)) %>%
  ungroup %>% 
  filter(!duplicated(paste0(size, decay))) %>% 
  arrange(desc(mean))
# # A tibble: 100 × 7
# Fold   size   decay   Auc  Sens  Spec  mean
# <fct> <int>   <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Fold1     6 0.01    0.955 0.893 0.922 0.951
# 2 Fold1     5 0.01    0.954 0.902 0.916 0.951
# 3 Fold1     3 0.0001  0.954 0.902 0.916 0.951
# 4 Fold1     3 0       0.952 0.902 0.919 0.951
# 5 Fold1     3 0.001   0.953 0.902 0.916 0.951
# 6 Fold1     6 0.0001  0.954 0.893 0.919 0.950
# 7 Fold1     3 0.00001 0.952 0.902 0.919 0.950
# 8 Fold1     6 0.00001 0.953 0.893 0.924 0.950
# 9 Fold1     4 0.001   0.952 0.902 0.919 0.950
# 10 Fold1     6 0.001   0.955 0.893 0.924 0.950
# # ℹ 90 more rows
# # ℹ Use `print(n = ...)` to see more rows

# Visulization----
summary.ml <- bind_rows(
  res$LR %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "LR"),
  res$NB %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "NB"),
  res$KNN %>% filter(para.k == 50) %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "KNN"),
  res$RF %>% filter(para.mtry == 1) %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "RF"),
  res$SVM %>% filter(para.cost == 1 & para.gamma == 100) %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "SVM"),
  res$XGBoost %>%
    group_by(Fold) %>% 
    filter(max_depth == 2 & eta == 0.2 & gamma == 0.01 & colsample_bytree == 0.5 & min_child_weight == 0.5) %>% 
    dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "XGBoost"),
  res$NN %>% filter(size == 6 & decay == 0.01) %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "NN"),
  res$DT %>% filter(para.cp == 0.0405) %>% dplyr::select(Fold, Auc, Sens, Spec) %>% mutate(ML = "NN")
) %>% 
  group_by(ML) %>% 
  mutate(mean = mean(Auc), sd = sd(Auc)) %>% 
  arrange(mean) 

summary.ml <- summary.ml %>% 
  mutate(ML = factor(ML, levels = unique(summary.ml$ML)))

p <- ggplot(summary.ml, aes(x = Auc, y = ML)) +
  geom_jitter(height = 0.3) +
  geom_point(data = summary.ml, aes(x = mean, y = ML), 
             size = 3, shape = 15, color = "#88cee6") +
  geom_errorbarh(data = summary.ml, aes(xmin = mean -sd, xmax = mean + sd, y = ML), 
                 width = .05, color = "#88cee6") +
  labs(y = "") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave(p, filename = "ML summary(Manually).pdf", width=8, height=8, units="in")