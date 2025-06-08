rm(list = ls())
gc()
save.image()

# load pkg----------------------------------------------------------------------
require(dplyr)
require(stringr)
require(sva)

# ComBat------------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/539507005
# 建立批次效应的模型，meta$group表示的是数据中除了有不同的批次，还有生物学上的差异。
# 在校正的时候要保留生物学上的差异，不能矫正过枉。

# ComBat_seq is an improved model from ComBat using negative binomial regression, which specifically targets RNA-Seq count data.
# sva::ComBat_Seq算法要求输入为count，对于以非Count形式给出的数据，由于无法转换回Count
# 参考35293897，可将RNA-seq数据转换为TPM，进而log2转换，以贴近芯片数据的分布形式

data.sva <- ComBat(mydata, batch = meta$GEO, 
                   mod = model.matrix(~as.factor(meta$group)))

pca_bf_sva <- prcomp(t(mydata), scale. = T, rank. = 3)

pca_af_sva <- prcomp(t(data.sva), scale. = T, rank. = 3)
# limma::removeBatchEffect------------------------------------------------------
