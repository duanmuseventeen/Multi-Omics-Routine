# 数据来源: 2024-Cell-38181739 Pozniak

# ----------------------------------------------------------
# 虚拟敲除
# ----------------------------------------------------------
# scTenifoldKnk----

library(scTenifoldKnk)

Malignant     = readRDS("Malignant_cells.rds")
Malignant_ds = subset(Malignant, downsample = 200)

count_mat <- SeuratObject::LayerData(Malignant_ds, assay = "RNA", layer = "counts")
count_mat <- as(count_mat, "dgCMatrix")

dim(count_mat)

ko_genes <- c(
  "ANTXR1",
  "CACNA1D",
  "CALCRL",
  "HGF",
  "LOXL1",
  "NGFR",
  "TRPV4"
)

gene_check <- data.frame(
  gene = ko_genes,
  in_matrix = ko_genes %in% rownames(count_mat)
)

gene_check$n_expr_cells <- NA_integer_
gene_check$pct_expr_cells <- NA_real_

genes_present <- intersect(ko_genes, rownames(count_mat))

gene_check$n_expr_cells[gene_check$gene %in% genes_present] <-
  Matrix::rowSums(count_mat[genes_present, , drop = FALSE] > 0)

gene_check$pct_expr_cells[gene_check$gene %in% genes_present] <-
  Matrix::rowMeans(count_mat[genes_present, , drop = FALSE] > 0) * 100

gene_check

hvg_use <- VariableFeatures(Malignant_ds)

network_genes <- unique(c(hvg_use, gene_check$gene))
network_genes <- intersect(network_genes, rownames(count_mat))

count_mat_use <- count_mat[network_genes, , drop = FALSE]

dim(count_mat_use)



outdir <- "scTenifoldKnk_results"
dir.create(outdir, showWarnings = FALSE)

ncores_use <- 16 # 问题1：linux系统节点中运行时，该参数无效

set.seed(1011)
scTenifoldKnk(
  countMatrix = count_mat_use,
  gKO = gene_check$gene[1],
  # transcriptomeWide = TRUE,
  
  qc = FALSE,
  
  nc_nNet = 10,
  nc_nCells = round(0.5 * ncol(count_mat_use)),
  nc_nComp = 3,
  nc_q = 0.9,
  
  td_K = 3,
  ma_nDim = 2,
  nCores = ncores_use
)

