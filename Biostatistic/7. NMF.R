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
# python------------------------------------------------------------------------
# https://mp.weixin.qq.com/s/YHgRZlun812jJDkf2fnB5g

conda activate python3.11
# python代码运行
# prepare步骤 
-c COUNTS, --counts COUNTS
[prepare] Input (cell x gene) counts matrix as .h5ad, .mtx,
df.npz, or tab delimited text file

cnmf prepare --output-dir ./res \ # 输出的文件夹为res
--name cNMF_res \ # 命名为cNMF_res
-c ./subdat.txt \ # 指定输入文件
-k 2 3 4 5 6 7 8 9 10 11 12 13 \ # 指定要尝试的因子数
--n-iter 100 --seed 123# 每个k值运行100次，并设定种子数

# factorize步骤
cnmf factorize --output-dir ./res \ # 指定输出目录，与prepare阶段一致，确保读取之前准备好的数据
--name cNMF_res \ # 与prepare保持一致，告诉程序使用哪一组预处理结果
--worker-index 0 --total-workers 1 # 当前工作进程的编号，编号从0开始。同时单核运行
# --total-workers必须为1，否则输出文件缺失

# combine步骤
cnmf combine --output-dir ./res --name cNMF_res # 执行合并步骤进行结果合并

# k_selection_图
cnmf k_selection_plot --output-dir ./res \ # 生成K值评估图，指定分析相同路径
--name cNMF_res # 

cnmf consensus --output-dir ./res \ # 执行共识分析，读取结果
--name cNMF_res \ # 与前面保持一致
--components 4 \ # 最终选定的因子数 k=4
--local-density-threshold 2 \ # 调整聚类时的密度阈值
--show-clustering # 输出可视化结果















