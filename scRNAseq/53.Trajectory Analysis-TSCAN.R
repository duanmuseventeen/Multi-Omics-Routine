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

# some cells cannot be assign a pseudotime
dim(lpsorder)
# [1] 5997    3
dim(fig)
# [1] 25414  7187
class(lpsorder)
# [1] "data.frame"

myorder <- lpsorder$Pseudotime
names(myorder) <- lpsorder$sample_name
# The tutorial is wrong
diffval <- difftest(procdata, myorder)
# take very long time...
diffval[diffval$qval < 0.05, ]

# Use singlegeneplot function to plot the expression value of a single gene against
# a given pseudotime. Notice that here orderonly should be set to FALSE.
STAT2expr <- log2(lpsdata["STAT2",]+1)
singlegeneplot(STAT2expr, TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))
# Comparing Different Pseudotemporal Ordering
# Not Run
# Where to use it?