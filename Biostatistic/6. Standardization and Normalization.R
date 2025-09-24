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
# Exploration 3-----------------------------------------------------------------
# outliers
# 异常值处理的重点是:
# 理解异常值的产生原因，因地制宜，采取对应的处理办法

# https://statorials.org/cn/box-cox-%E5%8F%98%E6%8D%A2%E5%88%B0-r/

library (MASS)

#create data
y=c(1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 6, 7, 8)
x=c(7, 7, 8, 3, 2, 4, 4, 6, 6, 7, 5, 3, 3, 5, 8)

#fit linear regression model
model <- lm(y~x)

#find optimal lambda for Box-Cox transformation 
bc <- boxcox(y ~ x)
(lambda <- bc$x[which.max(bc$y)])

[1] -0.4242424

#fit new linear regression model using the Box-Cox transformation
new_model <- lm(((y^lambda-1)/lambda) ~ x)

#define plotting area
op <- par(pty = "s", mfrow = c(1, 2))

#QQ plot for original model
qqnorm(model$residuals)
qqline(model$residuals)

#QQ plot for Box-Cox transformed model
qqnorm(new_model$residuals)
qqline(new_model$residuals)

#display both QQ plots
by(op)

# Exploration 4-----------------------------------------------------------------
# outliers
load("mymetabolomics.Rdata")

dat4heat <- dat.wide %>% 
  dplyr::select(-ID, -group, -PCA)

apply(dat4heat, 2, hist)
apply(log(dat4heat), 2, hist)

annot_col <- data.frame(
  PCA = dat.wide$PCA %>% factor,
  row.names = rownames(dat.wide)
)
# without treatment----
pheatmap::pheatmap(dat4heat,scale = "column", annotation_row = annot_col)

# scale df----
dat4heat %>% 
  scale %>% 
  pheatmap::pheatmap(scale = "none", annotation_row = annot_col)

# scale column----
tmp <- dat4heat
tmp[,1:ncol(tmp)] <- apply(tmp, 2, scale)

apply(tmp, 2, summary)

pheatmap::pheatmap(tmp, scale = "none", annotation_row = annot_col)

# IQR----
myout.IQR <- function(x){
  y <- x
  x[y < (quantile(y, .25) - 1.5 * IQR(y))] <- (quantile(y, .25) - 1.5 * IQR(y))
  x[y > (quantile(y, .75) + 1.5 * IQR(y))] <- (quantile(y, .75) + 1.5 * IQR(y))
  return(x)
}
tmp <- dat4heat
tmp[,1:ncol(tmp)] <- apply(tmp, 2, myout.IQR)

apply(tmp, 2, summary)

# pheatmap::pheatmap(tmp, scale = "column", annotation_row = annot_col)
pheatmap::pheatmap(scale(tmp), scale = "none", annotation_row = annot_col)
# 95%----
myout.95 <- function(x){
  y <- x
  x[y <= quantile(y, .0025)] <- quantile(y, .0025)
  x[y >= quantile(x, .9975)] <- quantile(y, .9975)
  return(x)
}
tmp <- dat4heat
tmp[,1:ncol(tmp)] <- apply(tmp, 2, myout.95)

apply(tmp, 2, summary)

# pheatmap::pheatmap(tmp, scale = "column", annotation_row = annot_col)
pheatmap::pheatmap(scale(tmp), scale = "none", annotation_row = annot_col)

# log----
dat4heat.log <- log(dat4heat + 1e-6)

pheatmap::pheatmap(dat4heat.log, scale = "none", annotation_row = annot_col)

dat4heat.log %>% 
  scale %>% 
  pheatmap::pheatmap(scale = "none", annotation_row = annot_col)

dat4heat.log %>% 
  pheatmap::pheatmap(scale = "column", annotation_row = annot_col)

tmp <- dat4heat.log
tmp[,1:ncol(tmp)] <- apply(tmp, 2, myout.IQR)
pheatmap::pheatmap(tmp, scale = "column", annotation_row = annot_col)











