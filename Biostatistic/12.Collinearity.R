# 方差膨胀因子（Variance Inflation Factor, VIF）
# 这是判断多重共线性最常用、最可靠的方法。VIF可以衡量一个自变量被其他所有自变量共同解释的程度。
# 计算方法： 对于模型中的每个自变量 X_j，VIF的计算方法是：
# 将 X_j 作为因变量，其他所有自变量作为自变量，运行一个线性回归模型。
# 计算这个回归模型的R2值。

# VIFj = 1/(1−Rj^2)

# 判断标准：
# VIF = 1：不存在共线性。
# VIF > 5：通常被认为是存在中等程度共线性的迹象。
# VIF > 10：通常被认为是存在严重共线性的迹象。
# VIF比简单的相关系数更强大，因为它能捕捉到一个变量与多个变量之间的关系。

setwd("C:/D/R project/Multi-Omics-Routine/Biostatistic/")

load("mymetabolomics.Rdata")

dat.wide <- dat.wide %>% dplyr::select(-group, -ID)

# lm----------------------------------------------------------------------------
fit <- lm(PCA ~ ., data = dat.wide)

car::vif(fit)
# Oleoylethanolamide         `Oleoyl-L carnitine`      `Palmitoyl-L carnitine` 
# 2.665268                    16.227739                     2.674828 
# `Stearoyl-L carnitine`        `Lauroyl-L carnitine` `cis-4-Decenoyl-L carnitine` 
# 14.590384                     4.952576                     6.358755 
# `Decanoyl-L carnitine`       `Hexanoyl-L carnitine`           `LysoPC(15:0/0:0)` 
# 12.344994                     4.893607                     8.253581 
# `LysoPC(17:0/0:0)`       `LysoPC(0:0/18:1(9Z))` 
# 11.241388                     3.890114 
# glm---------------------------------------------------------------------------
fit <- glm(PCA ~ ., data = dat.wide, family = binomial())

car::vif(fit)
# Oleoylethanolamide         `Oleoyl-L carnitine`      `Palmitoyl-L carnitine` 
# 2.979900                     5.205032                     2.631011 
# `Stearoyl-L carnitine`        `Lauroyl-L carnitine` `cis-4-Decenoyl-L carnitine` 
# 9.489220                    16.700980                     7.247455 
# `Decanoyl-L carnitine`       `Hexanoyl-L carnitine`           `LysoPC(15:0/0:0)` 
# 29.095244                    10.576281                    26.806176 
# `LysoPC(17:0/0:0)`       `LysoPC(0:0/18:1(9Z))` 
# 25.508579                     7.162776  
# Revision---------------------------------------------------------------------------
# 1. 保留一个最相关的变量
# 2. 结合或创造新变量
# 3. 使用正则化方法(e.g. lasso)
# 4. 主成分分析（PCA
fit <- glm(PCA ~ `Oleoyl-L carnitine` + `Stearoyl-L carnitine` + `LysoPC(17:0/0:0)`,
           data = dat.wide, family = binomial())
car::vif(fit)
# `Oleoyl-L carnitine` `Stearoyl-L carnitine`     `LysoPC(17:0/0:0)` 
# 1.303176               1.464711               1.220669





