# compare continuous and ordered categorical variables

t.test(...)

wilcox.test(...)

# multiple comparisons with pos-hoc---------------------------------------------
# one way ANOVA-----------------------------------------------------------------
dat <- readxl::read_excel("3. comparison.xlsx")

dat <- dat %>%
  mutate(group = paste0(group, "-", Trt))

# LSD ## 结果和Graphpad 9 一致 ## 结果和two-way ANOVA一致
require(DescTools)
tmp <- aov(value ~ group, dat)
res.lsd <- DescTools::PostHocTest(tmp, method = "lsd")

# HSD ## 结果和Graphpad 9 一致 ## 结果和two-way ANOVA一致
require(DescTools)
tmp <- aov(value ~ group, dat)
res.hsd <- DescTools::PostHocTest(tmp, method = "hsd")

# holm ## 结果和Graphpad 9 显著性一致，pvalue略有差别
res.holm_sidak <- pairwise.t.test(dat$value, dat$group, p.adjust.method = 'holm')

# newman-keuls ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group, dat)
res.nk <- DescTools::PostHocTest(tmp, method = "newmankeuls")

# Bonferroni ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group, dat)
res.bon <- DescTools::PostHocTest(tmp, method = "bonferroni")

# Dunnett's test 
# 结果与graphpad结果一致【仅适用于多组数据，均与其中一组比较的情况】
# https://search.r-project.org/CRAN/refmans/DescTools/html/DunnettTest.html
dat <- data.frame(
  group = c(rep("female-Ctrl",4),rep("female-MI",5),rep("female-OVX-MI",4)),
  value = c(10.34,4.33,4.87,14.45,
            9.27,18.74,18.18,18.29,11.61,
            14.72,10.19,15.58,21.72),
  stringsAsFactors = F
)

require(DescTools)
res.dunnett <- DescTools::DunnettTest(value ~ group, dat = dat, control = "female-Ctrl")

# two way ANOVA-----------------------------------------------------------------
dat <- readxl::read_excel("3. comparison.xlsx")

# LSD ## 未给出p value
require(agricolae)
tmp <- aov(value ~ group * Trt, dat)
res.lsd <- agricolae::LSD.test(tmp, trt=c("group","Trt"), p.adj = 'none')
res.snk <- agricolae::SNK.test(tmp, trt=c("group","Trt"),console = T)

# LSD ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group * Trt, dat)
res.lsd <- DescTools::PostHocTest(tmp, method = "lsd")

# HSD ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group * Trt, dat)
res.hsd <- DescTools::PostHocTest(tmp, method = "hsd")

# Tukey ## 结果和Graphpad 9 # (样本量不一致)无法计算
require(DescTools)
tmp <- aov(value ~ group * Trt, dat)
res.tukey <- TukeyHSD(tmp, "Group", ordered = TRUE)

# newman-keuls ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group * Trt, dat)
res.nk <- DescTools::PostHocTest(tmp, method = "newmankeuls")

# Bonferroni ## 结果和Graphpad 9 一致
require(DescTools)
tmp <- aov(value ~ group * Trt, dat)
res.bon <- DescTools::PostHocTest(tmp, method = "bonferroni")


# KW test-----------------------------------------------------------------
dat <- readxl::read_excel("3. comparison.xlsx")

dat <- dat %>%
  mutate(group = paste0(group, "-", Trt))

# Dunn
fit <- kruskal.test(value ~ group, data = dat)

library(dunn.test)
# fit1 <- dunn.test(dat$value, dat$group, method="bonferroni")

library(FSA)
fit2 <- dunnTest(value ~ group, data=dat, method="bonferroni") # selected

# Nemenyi test
require(DescTools)
fit <- NemenyiTest(value ~ group, data=dat) # selected





