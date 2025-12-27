library(rsvd)
library(scMetabolism)
library(Seurat)

# scMetabolism------------------------------------------------------------------
# https://github.com/wu-yc/scMetabolism

# 2. Quantify single-cell metabolism with Seurat (Recommended)
countexp.Seurat <- sc.metabolism.Seurat(
  obj = fig, # obj is a Seurat object containing the UMI count matrix.
  method = "AUCell", 
  imputation = F, # imputation allows users to choose whether impute their data before metabolism scoring.
  ncores = 4,
  metabolism.type = "KEGG")
# 错误于sc.metabolism.Seurat(obj = fig, method = "AUCell", imputation = F, :
# 没有名称为"counts"的插槽对于此对象类 "Assay5"

# To extract the metabolism score, just run metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score, where metabolism.matrix is the matrix.

# 4. Quantify single-cell metabolism WITHOUT Seurat (Not recommended)
# scMetabolism also supports quantifying metabolism independent of Seurat.
metabolism.matrix <- sc.metabolism(
  countexp = fig@assays$RNA$counts %>% as.data.frame,
  method = "AUCell",
  imputation = F,
  ncores = 4,
  metabolism.type = "KEGG")
# countexp is a data frame of UMI count matrix (col is cell ID, row is gene name).

fig@assays$METABOLISM$score <- metabolism.matrix
class(fig@reductions$umap@cell.embeddings)
# [1] "matrix" "array" 

countexp.Seurat <- fig
colnames(countexp.Seurat@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
# 3. Visualize
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)

input.pathway <- c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "drinking", norm = "y")

BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "drinking", ncol = 1)
