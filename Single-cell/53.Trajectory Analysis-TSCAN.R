library(TSCAN)
library(Seurat)

#============================== TSCAN ==========================================
# https://github.com/zji90/TSCAN
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://bioconductor.org/packages/devel/bioc/vignettes/TSCAN/inst/doc/TSCAN.pdf
# GUI---------------------------------------------------------------------------
TSCANui()
# Programming-------------------------------------------------------------------
procdata <- fig@assays$RNA$counts %>% 
  as.matrix %>% 
  preprocess(minexpr_value = 0, minexpr_percent = 0, cvcutoff = 0)
# The gene all express 0 will be removed.

# The function exprmclust will perform dimension reduction using principal component
# analysis and model-based clustering.
lpsmclust <- exprmclust(procdata)
# take very long time...
# Use plotmclust function to visualize the clustering based on the data after dimension
# reduction.
plotmclust(lpsmclust)
lpsorder <- TSCANorder(lpsmclust)
lpsorder

dim(lpsorder)
# [1] 5997    3
dim(fig)
# [1] 25414  7187
class(lpsorder)
# [1] "data.frame"

fig$id <- rownames(fig@meta.data)
fig.sub <- subset(fig, subset = id %in% lpsorder$sample_name)
meta <- fig@meta.data
lpsorder$id <- lpsorder$sample_name
meta <- meta %>% 
  left_join(lpsorder, by = "id")
meta <- as.data.frame(meta)
rownames(meta) <- meta$id
fig@meta.data <- as.data.frame(meta)
DimPlot(fig, group.by = "cell_type")
DimPlot(fig, group.by = "State")
FeaturePlot(fig, features = "Pseudotime")

myorder <- lpsorder$Pseudotime
names(myorder) <- lpsorder$sample_name
# The tutorial is wrong
diffval <- difftest(procdata, myorder)
diffval[diffval$FDR < 0.05, ]
#         pval           FDR
# PTGS2   0.000000e+00  0.000000e+00
# IL1R2   0.000000e+00  0.000000e+00
# HLA-B   0.000000e+00  0.000000e+00
# ACTB    0.000000e+00  0.000000e+00
# NAMPT  1.595832e-321 6.111399e-318
# S100A9 4.057465e-319 1.294872e-315

# Use singlegeneplot function to plot the expression value of a single gene against
# a given pseudotime. Notice that here orderonly should be set to FALSE.
PTGS2expr <- log1p(fig["PTGS2",])
singlegeneplot(PTGS2expr, TSCANorder(lpsmclust,orderonly=TRUE,flip=TRUE))
# 错误于exprdata$State: $ operator is invalid for atomic vectors
# Comparing Different Pseudotemporal Ordering
# Not Run
# Where to use it?