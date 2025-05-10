# https://www.jianshu.com/p/1cf65bd9aae7
# https://blog.csdn.net/jeffery0207/article/details/84348117

library(NMF)
library(dplyr)
library(Seurat)

# pmid: 40211000
# Deciphering the intratumoral transcriptional programs. To investigate
# potential intratumoral signatures of PCa, we performed NMF
# (cNMF package v.1.14 in Python) on malignant epithelia54. We first
# extracted and scaled the expression of the top 3,000 highly variable
# genes and replaced all negative expression values with 0. The top
# co-expressed gene modules in each sample were profiled, from which
# we extracted the 50 genes with the highest weight and used them to
# define intratumoral transcriptional programs. We obtained 3 to 20
# transcriptional programs in each sample and retained the ones with
# standard deviations larger than 0.2 and identified MPs shared among
# different samples through Jaccard index. Genes with expression percentages
# of more than 30% among all malignant epithelia were defined
# as hub genes for MPs.

# pmid: 37258682
# Defining a non-redundant set of robust NMF programs
# We carried out NMF for each sample separately, to generate programs
# that capture the heterogeneity within each sample. Negative values
# in each centred expression matrix were set to zero.

load("seurat.Rdata")

counts <- scobj.T.h@assays$RNA$counts[rownames(scobj.T.h@assays$RNA$scale.data),] %>%
  as.matrix
scale_data <- counts[,1:1000] %>% apply(1, scale) %>% t
colnames(scale_data) <- colnames(counts)[1:1000]
apply(scale_data, 1, range) %>% t

#  ...
# PLA2G2C             NaN      NaN
# C1QB         -0.1914854 7.467931
# LUZP1        -0.5198097 4.293243
# E2F2                NaN      NaN
# IFNLR1       -0.2159843 6.983492
# NCMAP               NaN      NaN
# CLIC4        -0.3545455 5.100000
# STMN1        -0.1421411 6.964912
# ZNF683       -0.1000000 9.900000
# TRNP1               NaN      NaN
# ...

scale_data <- scale_data[complete.cases(scale_data),]
scale_data[scale_data < 0] <- 0

# signature(x = "matrix", rank = "numeric", method = "NULL"): 
# Fits an NMF model using an appropriate algorithm when method is not supplied.
res_nmf <- nmf(scale_data, 5)
res_nmf
# <Object of class: NMFfit>
#   # Model:
#   <Object of class:NMFstd>
#   features: 1437 
# basis/rank: 5 
# samples: 100 
# # Details:
# algorithm:  brunet 
# seed:  random 
# RNG: 10403L, 10L, ..., -655188431L [8bf92d32e0bbbae79d3296b9443142f1]
# distance metric:  'KL' 
# residuals:  60166.17 
# Iterations: 1010 
# Timing:
#   用户 系统 流逝 
# 1.00 0.02 1.04 

extractFeatures(res_nmf)
# [[1]]
# <NA> 
#   NA 
# 
# [[2]]
# <NA> 
#   NA 
# 
# [[3]]
# <NA> 
#   NA 
# 
# [[4]]
# <NA> 
#   NA 
# 
# [[5]]
# <NA> 
#   NA 
# 
# attr(,"method")
# [1] "kim"

# 初次测试失败:
# 1. R-NMF包中NMF函数的参数和python中不同(以上两篇文献用的python包)
# 2. 使用降采样的数据运算，结果为NA，原因未知

# 之后尝试用Python测试









