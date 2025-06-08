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

# Example-----------------------------------------------------------------------
# This is a RNA sequencing data without known group of interest

# Note 1: if know length of genes, ComBat_seq are suggested for RNA seq data with counts as input
# Note 2: if know group of interest, it should be added into model.matrix
# Note 3: because the fpkm data corrected for batch effect here, we use ComBat fun
#         the data was log transformed to obtain a similar distribution to data from array

dat <- readxl::read_excel("All_gene_fpkm.xlsx")
dat <- dat %>% tibble::column_to_rownames("ID")

id <- colnames(dat)
batch <- id
batch[str_sub(id, 1, 1) %in% c("A","B","C")] <- 1
batch[str_sub(id, 1, 1) %in% c("D","E")] <- 2
batch <- as.factor(batch)

meta<- data.frame(
  # group = factor(group),
  row.names = id,
  stringsAsFactors = FALSE
)

dat.combat <- ComBat(
  dat = log2(dat + 1), 
  batch = batch, 
  mod = model.matrix(~1, data = meta), 
  par.prior=TRUE, 
  prior.plots=TRUE
)

# Found 4881 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
# Found2batches
# Adjusting for0covariate(s) or covariate level(s)
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data

dat.combat.exp <- exp(dat.combat) - 1

export::table2excel(dat.combat.exp, "expr_correted_batch.xlsx", add.rownames = TRUE)


