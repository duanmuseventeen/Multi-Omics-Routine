# https://www.jianshu.com/p/de25afe02a33
# https://www.cnblogs.com/yanzhi123/p/11712926.html
# pmid: 24766612
# pmid: 24942700
# pmid: 27694351

data <- mtcars

for (i in 1:ncol(data)) {
  tmp <- data[[i]]
  
  # mean normalization 
  data[[i]] <- tmp / mean(tmp)
  # median normalization 
  data[[i]] <- tmp / median(tmp)
  # max normalization 
  data[[i]] <- tmp / max(tmp)
  # global normalization (归一化)
  data[[i]] <- tmp / sum(tmp)
  # max-min normalization 
  data[[i]] <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
  # Z-score normaliztion
  data[[i]] <- (tmp - mean(tmp)) / sd(tmp)
  # IQR normaliztion
  data[[i]] <- (tmp - median(tmp)) / iqr(tmp)
}

# VSN
limma::normalizeVSN()

# Exploration 1-----------------------------------------------------------------
normalize <- function(x){
  return( (x - min(x)) / (max(x) - min(x)) )
}

a <- data$mpg %>% scale

b <- data$mpg %>% normalize

c <- data$mpg %>% normalize %>% scale

d <- data$mpg %>% scale %>% normalize

# 结论:
# a 等于 c，b 等于 d，因此连续进行scale和normalize是无意义的
# Exploration 2-----------------------------------------------------------------
err1 <- runif(32)
err3 <- err2 <- err1
err2[32] <- 1e6
err3[29:32] <- c(1e6, 0.74e5, 3.14e7, 0.22e3)

data.frame(
  a = data$mpg %>% scale,
  b = (data$mpg * 10) %>% scale,
  c = (data$mpg * 10 - err1) %>% scale,
  d = (data$mpg * 10 + err2) %>% scale,
  e = (data$mpg * 10 - err3) %>% scale
)

data.frame(
  a = data$mpg %>% mean,
  b = (data$mpg * 10) %>% mean,
  c = (data$mpg * 10 - err1) %>% mean,
  d = (data$mpg * 10 + err2) %>% mean,
  e = (data$mpg * 10 - err3) %>% mean
)

# 结论:
# scale [Z-score] 可以较好的消除量纲的影响
# scale [Z-score] 对异常值敏感，处理前一定要先处理异常值




