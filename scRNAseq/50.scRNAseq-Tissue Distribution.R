# Ref:
# https://github.com/neurorestore/Augur
# https://github.com/MarioniLab/miloR

# https://zhuanlan.zhihu.com/p/715794329
# https://blog.csdn.net/zfyyzhys/article/details/147021127

require(Seurat)
require(ggplot2)
require(dplyr)
require(Startrac)
require(miloR)
require(Augur)

load("Seurat.Rdata") # seurat.obj
# Ro/e =========================================================================
# tissue migration and state transition. 
# We present STARTRAC as a framework, defined by four indices, to analyse different 
# aspects of T cells based on paired single-cell transcriptomes and TCR sequences. 
# The first index, STARTRAC-dist, uses the ratio of observed over expected cell numbers 
# in tissues to measure the enrichment of T cell clusters across different tissues. 
# Given a contingency table of T cell clusters by tissues, we first apply chi-squared 
# test to evaluate whether the distribution of T cell clusters across tissues 
# significantly deviates from random expectations. We then calculate the STARTRAC-dist 
# index for each combination of T cell clusters and tissues according the following formula:

# r_STARTRAC_dist = R_o/e = observed / expected

# in which R_o/e is the ratio of observed cell number over the expected cell number 
# of a given combination of T cell cluster and tissue. The expected cell number 
# for each combination of T cell clusters and tissues are obtained from the 
# chi-squared test. Different from the chi-squared values, which are defined as 
# (observed - expected)^2 / expected and can only indicate the divergence of 
# observations from random expectations, r_STARTRAC_dist defined by R_o/e can 
# indicate whether cells of a certain T cell cluster are enriched or depleted in
# a specific tissue. For example, if R_o/e > 1, it suggests that cells of the 
# given T cell cluster are more frequently observed than random expectations in 
# the specific tissue, that is, enriched. If R_o/e < 1, it suggests that cells of 
# the given T cell cluster are observed with less frequency than random expectations 
# in the specific tissue, that is, depleted. By calculating the STARTRAC-dist indices 
# via R_o/e, we can quantify the tissue preference of T cell clusters efficiently.

data <- seurat.obj@meta.data
data$majorCluster = data$cell_type
data$patient = data$sample
data$loc = data$condition
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe

# col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), 
#                       c("blue", "white", "red"))
# Heatmap(as.matrix(Roe),
#         show_heatmap_legend = TRUE, 
#         cluster_rows = TRUE, 
#         cluster_columns = TRUE,
#         row_names_side = 'right', 
#         show_column_names = TRUE,
#         show_row_names = TRUE,
#         # col = col_fun,
#         row_names_gp = gpar(fontsize = 10),
#         column_names_gp = gpar(fontsize = 10),
#         heatmap_legend_param = list(
#           title = "Ro/e Index",  # 自定义图注名称
#           at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
#           labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
#         ),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
#         }
# )
# # cell_fun = function(j, i, x, y, width, height, fill)
# # j 和 i：分别表示单元格的列和行索引。
# # x 和 y：是单元格的中心位置坐标，用于确定文本的位置。
# # width 和 height：表示单元格的宽度和高度，可以用来调整文本位置或大小。
# # fill：表示单元格的背景颜色。
# MiloR ========================================================================
# https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html

# ...
# From Seurat object
# The Seurat package includes a converter to SingleCellExperiment.

seurat.obj_sce <- as.SingleCellExperiment(seurat.obj)
seurat.obj_milo <- Milo(seurat.obj_sce)
  
# 1. Build a graph and neighbourhoods.------------------------------------------
# We need to add the KNN graph to the Milo object. This is stored in the graph slot, 
# in igraph format. The miloR package includes functionality to build and store
# the graph from the PCA dimensions stored in the reducedDim slot.

# d	
# The number of dimensions to use if the input is a matrix of cells X reduced dimensions. If this is provided, transposed should also be set=TRUE.
milo.obj <- buildGraph(seurat.obj_milo, k=20, d=20)
# 2. Defining representative neighbourhoods-------------------------------------
milo.obj <- makeNhoods(milo.obj, k=20, d=20, refined=TRUE, prop=0.2)
plotNhoodSizeHist(milo.obj)
# 3. Counting cells in neighbourhoods-------------------------------------------
require(SingleCellExperiment)
milo.obj <- countCells(milo.obj, 
                       meta.data = data.frame(colData(milo.obj)), 
                       samples="sample")
head(nhoodCounts(milo.obj))
# 6 x 10 sparse Matrix of class "dgCMatrix"
# [[ suppressing 10 column names ‘T01’, ‘T02’, ‘T03’ ... ]]
# 1 1 2 11  . 2  .  2 23 11  .
# 2 . 2  6  . 6 10 12  7 14  1
# 3 . .  2 20 5  1  .  9  4  .
# 4 . .  5  . 1  7  6  5 16  .
# 5 3 1  3  . 3  3 15  5 13 13
# 6 . .  7  . 7  .  .  7  7  .
# 4. Differential abundance testing---------------------------------------------
design <- data.frame(colData(milo.obj))[,c("sample", "condition")]
design <- distinct(design)
rownames(design) <- design$sample
## Reorder rownames to match columns of nhoodCounts(milo)
design <- design[colnames(nhoodCounts(milo.obj)), , drop=FALSE]

design

milo.obj <- calcNhoodDistance(milo.obj, d=20)
## 'as(<dgTMatrix>, "dgCMatrix")' is deprecated.
## Use 'as(., "CsparseMatrix")' instead.
## See help("Deprecated") and help("Matrix-deprecated").

# rownames(design) <- design$sample
da_results <- testNhoods(milo.obj, design = ~ condition, design.df = design)

# This calculates a Fold-change and corrected P-value for each neighbourhood, 
# which indicates whether there is significant differential abundance between 
# conditions.
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 
#              logFC    logCPM            F    PValue       FDR Nhood SpatialFDR
# 9    -0.0007365177  9.948485 1.269405e-06 0.9991011 0.9991011     9  0.9991011
# 506  -0.0063837016 10.145997 2.431591e-05 0.9960656 0.9974539   506  0.9975369
# 605  -0.0109651231  9.740627 3.761177e-05 0.9951068 0.9974539   605  0.9975369
# 886   0.0071344580 10.313033 1.604091e-05 0.9968045 0.9974986   886  0.9975369
# 1223 -0.0098064441 10.454771 2.886528e-05 0.9957134 0.9974539  1223  0.9975369
# 168   0.0124211078  9.947118 7.345623e-05 0.9931618 0.9969152   168  0.9971294
# 5. Visualize neighbourhoods displaying DA-------------------------------------
milo.obj <- buildNhoodGraph(milo.obj)

scater::plotUMAP(milo.obj) + plotNhoodGraphDA(milo.obj, da_results, alpha=0.05) +
  patchwork::plot_layout(guides="collect")
# Augur ========================================================================
seurat.obj@meta.data$labels <- seurat.obj@meta.data$condition

augur <- calculate_auc(seurat.obj, cell_type_col = "cell_type", label_col = "condition")

# odds ratio ===================================================================
pmid: 37696831



augur$AUC

