## ML Model---------------------------------------------------------------------
# install.packages("sampling") # sample
# install.packages("xgboost") # xgboost
# install.packages("e1071") # NaiveBayes, SVM
# install.packages("kernlab") # SVM
# install.packages("klaR") # NaiveBayes
# install.packages("randomForest") # RandomForest
# install.packages("ranger") # RandomForest
# install.packages("neuralnet") # NNET
# install.packages("class") # KNN
# install.packages("FNN") # KNN
# install.packages("kknn") # KNN
# install.packages("rpart") # DT
# install.packages("caret") # cross validation & tune
suppressMessages(require(dplyr))
suppressMessages(require(ggplot2))
suppressMessages(require(stringr))

# Step 1: 按照7:3的比例，将数据拆分为内部训练和内部测试集
if(all(meta_GSE211692$title == colnames(expr_GSE211692))){
    grp_GSE211692 <- ifelse(meta_GSE211692$`disease state:ch1` == "no cancer", 0, 1)
    group <- grp_GSE211692
    tmp <- expr_GSE211692[rownames(expr_GSE211692) %in% up.core, ] %>% t
    dat.raw <- tmp %>% as.data.frame %>% 
      mutate(group = factor(group))
    colnames(dat.raw)[1:2] <- c('hsa_miR_1247_3p', 'hsa_miR_744_5p')
  } 
  
table(dat.raw$group)
# 0    1 
# 5643 1596

set.seed(seed)
# sample method 1----------------------------------------------------------------------
# ? Poisson sampling (poisson), systematic sampling (systematic)
# dat.sam <- strata(dat.raw %>% arrange(group), 'group',
#        size = round(table(dat.raw$group)*0.5), method = "srswor")
# sample method 2----------------------------------------------------------------------
# strat sampling
sele <- c(sample(rownames(dat.raw)[dat.raw$group == 0], size = 2822, replace = FALSE),
          sample(rownames(dat.raw)[dat.raw$group == 1], size = 798, replace = FALSE))

dat.find <- dat.raw[sele,]
dat.vali <- dat.raw[!rownames(dat.raw) %in% sele, ]

# k-fold sampling
piece <- round(nrow(dat.raw) / 10)
dataList <- list()
for (i in 1:10) {
  if((piece * i) >= nrow(dat.raw)){
    dataList[[i]] <- dat.raw[(piece * (i - 1) + 1):nrow(dat.raw),]
  }else{
    dataList[[i]] <- dat.raw[(piece * (i - 1) + 1):(piece * i),]
  }
}
# sample method 3----------------------------------------------------------------------
in_training <- createDataPartition(dat.raw$group, p = .80, list = FALSE)
training_data <- dat.raw[in_training, ]
testing_data  <- dat.raw[-in_training, ]

# Step 2: 在内部训练集中，使用嵌套交叉验证（蒙特卡洛或10折），评估不同模型的预测能力----
# https://topepo.github.io/caret/available-models.html

myML <- function(dat.find = dat.find, dat.vali = dat.vali, seed = 1011, seed4caret = c(1:10, 1011)){
  suppressMessages(require(e1071))
  suppressMessages(require(randomForest))
  suppressMessages(require(neuralnet))
  suppressMessages(require(class))
  suppressMessages(require(rpart))
  suppressMessages(require(ranger))
  suppressMessages(require(caret))
  suppressMessages(require(kernalshap))
  # suppressMessages(require(FNN))
  
  set.seed(seed)
  
  # 0. LR-------------------------------------------------------------------------
  # 10-fold sampling
  # piece <- round(nrow(dat.raw) / 10)
  # dataList <- list()
  # dat.left <- dat.raw
  # for (i in 1:10) {
  #   if(piece < nrow(dat.left)){
  #     mysam <- sample(rownames(dat.left), size = piece, replace = FALSE)
  #     dataList[[i]] <- dat.left[rownames(dat.left) %in% mysam, ]
  #     dat.left <- dat.left[!(rownames(dat.left) %in% mysam),]
  #   }else{
  #     dataList[[i]] <- dat.left
  #   }
  # }
  # for (i in 1:10) {
  #   dat.test.lr <- dataList[[i]]
  #   dat.find.lr <- dat.raw[!(rownames(dat.raw) %in% rownames(dataList[[i]])),]
  #   
  #   fit.lr <- glm(group ~ ., data = as.data.frame(dat.find.lr), family = binomial())
  #   modelroc.lr <- roc(dat.test.lr$group, predict(fit.lr, newdata = dat.test.lr[,1:2], type='response'))
  #   ci(modelroc.lr)
  # }

  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         preProcOptions = list(method = c("center", "scale")),
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(1:10, seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  lr.train = train(group ~ ., data = dat.find, method = "glm", metric = "ROC",
                   family = binomial(),
                   trControl = control)
  # lr.train
  # lr.train$resample
  # 1. {scale-insens} NB[Tune]-------------------------------------------------------------------------
  # require(e1071)
  # fit.nb <- naiveBayes(group ~ ., data = dat.find, laplace=0)
  # modelroc.nb <- roc(dat.vali$group, predict(fit.nb, dat.vali) %>% as.numeric)
  # ci(modelroc.nb)

  # require(klaR)
  # fit.NB <- NaiveBayes(group ~ ., data = dat.find)
  # # plot(fit.NB)
  # modelroc.NB <- roc(dat.vali$group, predict(fit.NB, dat.vali)$class %>% as.numeric)
  # ci(modelroc.NB)
  
  grid = expand.grid(usekernel = c(FALSE),
                     fL = c(0, 0.5, 1),
                     adjust = 1)
  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(lapply(as.list(rep(100,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                  seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  nb.train = train(group ~ ., data = dat.find, method = "nb", metric = "ROC",
                   trControl = control, tuneGrid = grid)
  # nb.train
  # nb.train$resample
  # 2. {tree} XGBoost[Tune]--------------------------------------------------------------------
  # https://xgboost.readthedocs.io/en/release_1.5.0/R-package/xgboostPresentation.html
  # Extreme Gradient Boosting, which is an efficient implementation of the gradient boosting framework from Chen & Guestrin (2016) <doi:10.1145/2939672.2939785>. This package is its R interface. The package includes efficient linear model solver and tree learning algorithms. The package can automatically do parallel computation on a single machine which could be more than 10 times faster than existing gradient boosting packages. It supports various objective functions, including regression, classification and ranking. The package is made to be extensible, so that users are also allowed to define their own objectives easily.
  
  # bst <- xgboost(data = dat.find[,1:2] %>% as.matrix, 
  #                label = dat.find$group %>% as.character %>% as.numeric, 
  #                max.depth = 3, # 通常 3 ~ 6
  #                eta = 0.1, # 通常 0.01 ~ 0.3
  #                nrounds = 10, colsample_bytree=1,
  #                nthread = 4, 
  #                objective = "binary:logistic" # 二分类问题
  # )
  # # Tune
  # 
  # modelroc.bst <- roc(dat.vali$group, predict(bst, dat.vali[,1:2] %>% as.matrix))
  # ci(modelroc.bst)
  
  grid = expand.grid(nrounds = 1000,
                     max_depth = c(3:6),
                     eta = c(0.1, 0.2, 0.3),
                     gamma = c(0.01, 0.1, 0.2, 0.5, 1.0),
                     colsample_bytree = 1,
                     min_child_weight = c(0.5, 1, 1.5),
                     subsample = 1)
  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(lapply(as.list(rep(1000,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                  seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  xgb.train = train(group ~ ., data = dat.find, method = "xgbTree", metric = "ROC",
                    # objective = "binary:logistic", 
                    nthread = 4,
                    trControl = control, tuneGrid = grid)
  # xgb.train
  # xgb.train$resample
  # 3. {tree} RF [Tune]-------------------------------------------------------------------------
  # fit.rf <- randomForest(group ~ ., data=dat.find, ntree=500, mtry=2,
  #                        importance=FALSE, proximity=FALSE)
  # modelroc.rf <- roc(dat.vali$group, predict(fit.rf, dat.vali) %>% as.numeric)
  # 
  # # fit.RF <- ranger::ranger(group ~ ., data = dat.find, num.trees = 1000, mtry = 1, importance = "none", probability = TRUE)
  # # modelroc.RF <- roc(dat.vali$group, predict(fit.RF, dat.vali[,1:2])$predictions %>% as.numeric)
  
  grid = expand.grid(mtry = seq(1,2, by = 1))
  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(lapply(as.list(rep(100,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                  seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  rf.train = train(group ~ ., data = dat.find, method = "rf", metric = "ROC",
                   ntree=500, importance=FALSE, proximity=FALSE,
                   trControl = control, tuneGrid = grid)
  # rf.train
  # rf.train$resample
  # 4. {scale-sens} SVM------------------------------------------------------------------------
  # fit.svm <- e1071::svm(group ~ ., data=dat.find, kernel = "radial")
  # tuneResult <- tune(svm, group ~ ., data=dat.find, kernel = "radial",
  #                    ranges = list(cost = c(0.1, 1, 10), gamma = c(0.1, 1, 10)))
  # bestModel <- tuneResult$best.model
  # modelroc.svm <- roc(dat.vali$group, predict(bestModel, dat.vali[,1:2]) %>% as.numeric)
  
  # fit.ksvm <- kernlab::ksvm(group ~ ., data = dat.find, type = "C-svc", 
  #                           kernel="rbfdot", # Radial Basis kernel "Gaussian"
  #                           kpar=list(sigma=0.05), #
  #                           C = 1, # 
  #                           cross = 3, 
  #                           na.action = na.omit, scaled = FALSE)
  # modelroc.ksvm <- roc(dat.vali$group, predict(fit.ksvm, dat.vali[,1:2]) %>% as.numeric)
  # modelroc.ksvm
  
  grid = expand.grid(sigma = c(0.01, 0.1, 1, 10, 100),
                     C = c(0.01, 0.1, 1, 10, 100))
  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(lapply(as.list(rep(100,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                  seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  svm.train = train(group ~ ., data = dat.find,
                    method = "svmRadial",  metric = "ROC", 
                    na.action = na.omit, scaled = FALSE,
                    trControl = control, tuneGrid = grid)
  # svm.train
  # svm.train$resample
  # 5. {scale-sens} KNN [Tune]------------------------------------------------------------------------
    # knn_pred <- knn(train = dat.find[,1:2], test = dat.vali[,1:2], 
    #                 cl = dat.find$group, k = k)
    # modelroc.knn <- roc(dat.vali$group, knn_pred %>% as.numeric)
    # ci(modelroc.knn)
    
    # kknn_pred <- kknn::kknn(formula = group ~ ., dat.find, dat.vali, na.action = na.omit(),
    #                         k = 2, scale=TRUE, kernel = "optimal",
    #                         distance = 2, # 2 代表欧几里得距离
    #                         ykernel = NULL,
    #                         contrasts = c('unordered' = "contr.dummy", ordered = "contr.ordinal"))
    # modelroc.kknn <- roc(dat.vali$group, kknn_pred$fitted.values %>% as.numeric)
    # ci(modelroc.knn)
    
    # fnn_pred <- FNN::knn(train = dat.find[,1:2], test = dat.vali[,1:2], cl = dat.find$group, k = 2)
    # modelroc.fnn <- roc(dat.vali$group, fnn_pred %>% as.numeric)
    # ci(modelroc.knn)
  
    grid = expand.grid(.k = seq(5, 50, by = 5))
    set.seed(seed)
    control = trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           seed = c(lapply(as.list(rep(1000,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                    seed),
                           returnResamp = "all",
                           allowParallel = FALSE)
    knn.train = train(group ~ ., data = dat.find, method = "knn", metric = "ROC",
                      trControl = control, tuneGrid = grid)
    # knn.train
    # knn.train$results
    
    # knn.train$resample
    # knn.train$resample %>% group_by(k) %>% 
    #   mutate(Accuracy_SD = sd(Accuracy),
    #          Kappa_SD = sd(Kappa)) %>% 
    #   mutate(Accuracy = mean(Accuracy),
    #          Kappa = mean(Kappa)) %>% 
    #   distinct(k, .keep_all = TRUE)
    # 
    # knn.train$resample %>% group_by(Resample) %>% 
    #   filter(Accuracy == max(Accuracy))
    #   distinct(Resample, .keep_all = TRUE)
    # 
    # knn.train$control$index # 抽样结果
  # 6. {scale-sens} NNET-----------------------------------------------------------------------
    # single layer
    # fit.nn <- nnet::nnet(group ~ ., data = dat.find, size = 10, rang = 0.5,
    #                decay = 5e-4, maxit = 500, linout = FALSE)
    # modelroc.nn <- roc(dat.vali$group, predict(fit.nn, dat.vali[,1:2])[,1] %>% as.numeric)
    # ci(modelroc.nn)
    # NeuralNetTools::plotnet(fit.nn)
    
    # mulitiple layers
    # fit.neun <- neuralnet::neuralnet(
    #   group == 1 ~ ., dat.find, rep = 1, hidden = c(3,3),
    #   linear.output = FALSE)
    # modelroc.neun <- roc(dat.vali$group, predict(fit.neun, dat.vali[,1:2])[,1] %>% as.numeric)
    # ci(modelroc.neun)
    # plot(fit.neun)
    # NeuralNetTools::plotnet(fit.neun)
    
    grid = expand.grid(size = c(1:20),
                       decay = c(0, 1e-5, 1e-4, 1e-3, 1e-2))
    set.seed(seed)
    control = trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           seed = c(lapply(as.list(rep(1000,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                    seed),
                           returnResamp = "all",
                           allowParallel = FALSE)
    nn.train = train(group ~ ., data = dat.find %>% mutate(group = factor(group, labels = c("NC","CRC"))),
                     method = "nnet", metric = "ROC",
                     rang = 0.5, maxit = 1000, linout = FALSE,
                     trControl = control, tuneGrid = grid)
    # nn.train
    # nn.train$resample
  # 7. {tree} DT [Tune]-------------------------------------------------------------------------
  # fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  # modelroc.dt <- roc(dat.vali$group, predict(fit.dt, dat.vali[,1:2], type="class") %>% as.numeric)
  # ci(modelroc.dt)
  
  # fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  # fit.dt.tune <- rpart(group~., data=dat.find,
  #                      control=rpart.control(minsplit=30, minbucket = 10, cp=0.01, maxdepth = 30))
  # modelroc.dt <- roc(dat.vali$group, predict(fit.dt.tune, dat.vali[,1:2], type="class") %>% as.numeric)
  # ci(modelroc.dt)
  
  # fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  # # printcp(fit.dt)
  # # plotcp(fit.dt)
  # # fit.dt$cptable
  # cp <- fit.dt$cptable[,"xerror"]
  # ncp <- length(cp)
  # if(length(which((cp[1:ncp-1] - cp[2:ncp]) <= 0)) == 0){
  #   best_cp <- fit.dt$cptable[which.min(fit.dt$cptable[,"xerror"]), "CP"]
  # }else{
  #   best_cp <- fit.dt$cptable[which((cp[1:ncp-1] - cp[2:ncp]) <= 0)[1], "CP"]
  # }
  # fit.dt.prune <- prune(fit.dt, cp = best_cp)
  # modelroc.dt <- roc(dat.vali$group, predict(fit.dt.prune, dat.vali[,1:2], type="class") %>% as.numeric)
  # ci(modelroc.dt)
    
  grid = expand.grid(cp = seq(0.0001, .05, by = 0.0001))
  set.seed(seed)
  control = trainControl(method = "cv", 
                         number = 10, 
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         seed = c(lapply(as.list(rep(1000,10)), function(x){sample(c(1:x), nrow(grid), replace = FALSE)}),
                                  seed),
                         returnResamp = "all",
                         allowParallel = FALSE)
  dt.train = train(group ~ ., data = dat.find, method = "rpart", metric = "ROC",
                   trControl = control, tuneGrid = grid)
  # dt.train
  # dt.train$resample
  # 8. RES-------------------------------------------------------------------------
  res <- list(LR = lr.train, 
              # NB.e1071 = modelroc.nb$auc,
              # NB.klaR = modelroc.NB$auc,
              NB = nb.train,
              XGBoost = xgb.train,
              KNN = knn.train,
              RF = rf.train,
              SVM = svm.train,
              KNN = knn.train,
              NN = nn.train,
              DT = dt.train)
  save(res, file = "Flow/ML(scale).Rdata")
  
  return(list(LR = lr.train, 
              # NB.e1071 = modelroc.nb$auc,
              # NB.klaR = modelroc.NB$auc,
              NB = nb.train,
              XGBoost = xgb.train,
              KNN = knn.train,
              RF = rf.train,
              SVM = svm.train,
              KNN = knn.train,
              NN = nn.train
              DT = dt.train))
}


