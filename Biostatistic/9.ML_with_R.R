# install.packages("sampling") # sample
# install.packages("xgboost") # xgboost
# install.packages("e1071") # NaiveBayes
# install.packages("klaR") # NaiveBayes
# install.packages("randomForest") # RandomForest
# install.packages("ranger") # RandomForest
# install.packages("neuralnet") # NNET
# install.packages("class") # KNN
# install.packages("FNN") # KNN
# install.packages("kknn") # KNN
# install.packages("rpart") # DT
suppressMessages(require(dplyr))
suppressMessages(require(ggplot2))
suppressMessages(require(stringr))

table(dat.raw$group)

mymachiinelearning <- function(dat, size.finding.nc = 2822, size.finding.dis = 798, seed){
  suppressMessages(require(e1071))
  suppressMessages(require(randomForest))
  suppressMessages(require(neuralnet))
  suppressMessages(require(class))
  suppressMessages(require(rpart))
  suppressMessages(require(ranger))
  # suppressMessages(require(FNN))
  
  set.seed(seed)
  # ? Poisson sampling (poisson), systematic sampling (systematic)
  # dat.sam <- strata(dat.raw %>% arrange(group), 'group',
  #        size = round(table(dat.raw$group)*0.5), method = "srswor")
  
  sele <- c(sample(rownames(dat.raw)[dat.raw$group == 0], size = 2822, replace = FALSE),
            sample(rownames(dat.raw)[dat.raw$group == 1], size = 798, replace = FALSE))
  
  dat.find <- dat.raw[sele,]
  dat.vali <- dat.raw[!rownames(dat.raw) %in% sele, ]
  
  set.seed(1011)
  # 0. LR-------------------------------------------------------------------------
  fit.lr <- glm(group ~ ., data = as.data.frame(dat.find), family = binomial())
  modelroc.lr <- roc(dat.vali$group, predict(fit.lr, newdata = dat.vali, type='response'))
  ci(modelroc.lr)
  # 1. XGBoost--------------------------------------------------------------------
  # bst <- xgboost(data = dat.find, label = dat.find$group, max_depth = 2, eta = 1,
  #                nrounds = 2, objective = "binary:logistic", nthread = 2)
  # pred <- predict(bst, test$data)
  # 2. RF-------------------------------------------------------------------------
  fit.rf <- randomForest(dat.find$group ~ ., data=dat.find, ntree=500, mtry=2,
                         importance=FALSE, proximity=FALSE)
  modelroc.rf <- roc(dat.vali$group, predict(fit.rf, dat.vali) %>% as.numeric)
  
  # fit.RF <- ranger::ranger(dat.find$group ~ ., data = dat.find, num.trees = 500, importance = "none")
  # modelroc.RF <- roc(dat.vali$group, predict(fit.RF, dat.vali[,1:2])$predictions %>% as.numeric)
  # 3. SVM------------------------------------------------------------------------
  fit.svm <- e1071::svm(dat.find$group ~ ., data=dat.find, kernel = "sigmoid")
  tuneResult <- tune(svm, dat.find$group ~ ., data=dat.find, kernel = "sigmoid",
                     ranges = list(cost = c(0.1, 1, 10), gamma = c(0.1, 1, 10)))
  bestModel <- tuneResult$best.model
  
  predicted <- predict(bestModel, testData)
  # 4. KNN------------------------------------------------------------------------
  knn_pred <- knn(train = dat.find[,1:2], test = dat.vali[,1:2], cl = dat.find$group, k = 2)
  modelroc.knn <- roc(dat.vali$group, knn_pred %>% as.numeric)
  ci(modelroc.knn)
  
  # kknn_pred <- kknn::kknn(formula = group ~ ., dat.find, dat.vali, na.action = na.omit(), 
  #      k = 2, scale=TRUE, kernel = "optimal", 
  #      distance = 2, # 2 代表欧几里得距离 
  #      ykernel = NULL, 
  #      contrasts = c('unordered' = "contr.dummy", ordered = "contr.ordinal"))
  # modelroc.kknn <- roc(dat.vali$group, kknn_pred$fitted.values %>% as.numeric)
  # ci(modelroc.knn)
  
  # fnn_pred <- FNN::knn(train = dat.find[,1:2], test = dat.vali[,1:2], cl = dat.find$group, k = 2)
  # modelroc.fnn <- roc(dat.vali$group, fnn_pred %>% as.numeric)
  # ci(modelroc.knn)
  # 5. NNET-----------------------------------------------------------------------
  
  # 6. NB-------------------------------------------------------------------------
  require(e1071)
  fit.nb <- naiveBayes(group ~ ., data = dat.find, laplace=0)
  modelroc.nb <- roc(dat.vali$group, predict(fit.nb, dat.vali) %>% as.numeric)
  ci(modelroc.nb)
  
  require(klaR)
  fit.NB <- NaiveBayes(group ~ ., data = dat.find)
  modelroc.NB <- roc(dat.vali$group, predict(fit.NB, dat.vali)$class %>% as.numeric)
  ci(modelroc.NB)
  # 7. DT-------------------------------------------------------------------------
  # fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  # modelroc.dt <- roc(dat.vali$group, predict(fit.dt, dat.vali[,1:2], type="class") %>% as.numeric)
  # ci(modelroc.dt)
  
  fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  fit.dt.tune <- rpart(group~., data=dat.find,
                       control=rpart.control(minsplit=30, minbucket = 10, cp=0.01, maxdepth = 30))
  modelroc.dt <- roc(dat.vali$group, predict(fit.dt.tune, dat.vali[,1:2], type="class") %>% as.numeric)
  ci(modelroc.dt)
  
  # fit.dt <- rpart::rpart(group~., data=dat.find, method = "class")
  # # printcp(fit.dt)
  # # plotcp(fit.dt)
  # # fit.dt$cptable
  # best_cp <- fit.dt$cptable[which.min(fit.dt$cptable[,"xerror"]), "CP"]
  # fit.dt.prune <- prune(fit.dt, cp = best_cp)
  # modelroc.dt <- roc(dat.vali$group, predict(fit.dt.prune, dat.vali[,1:2], type="class") %>% as.numeric)
  # ci(modelroc.dt)
  # 8. RES-------------------------------------------------------------------------
  return(list(LR = modelroc.lr, 
              RF = modelroc.rf, 
              NB.e1071 = modelroc.nb,
              NB.klaR = modelroc.NB,
              DT = modelroc.dt))
}
