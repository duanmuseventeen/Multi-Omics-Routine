# Formula ---------------------------------------------------------------------
p_value <- (T_perm >= T_obs) / n 
or 
p_value <- ((T_perm >= T_obs) + 1) / (n + 1) 

# example 1 -------------------------------------------------------------------
B <- 1000
obs_auc <- .8

permutation.aucs <- lapply(c(1:B), function(x){
  set.seed(x)
  print(x)

  group <- sample(meta$group)
  fit <- cv.glmnet(as.matrix(t(data)), group, family = binomial(link = 'logit'), type.measure = 'auc')
  prob <- predict(fit, newx = as.matrix(t(data)), s = "lambda.min", type = "response")
  modelroc <- roc(group, unlist(prob[,1])) # label used here is not true label
  return(modelroc)
  })

aucs <- lapply(permutation.aucs, function(x) x$auc %>% as.numeric()) %>% unlist
p_value <- (sum(aucs >= obs_auc) + 1) / (B + 1)
