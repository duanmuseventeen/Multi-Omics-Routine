# Ref:
# https://github.com/digitalcytometry/cytotrace2
# https://mp.weixin.qq.com/s/KGSoRx3klmliKPVL7ml28Q

require(CytoTRACE2)
require(dplyr)
require(ggplot2)
require(Seurat)

load("scobj.Rdata")

# CytoTRACE 2 requires a single-cell RNA-sequencing gene expression object as input. This can include either raw counts or CPM/TPM normalized counts, and should not be log-transformed. This input can be provided in any of the following formats:
#   
# A data table. This object of class data.frame or of another class which allows storing row and column names. This should have genes as rows and cells as columns, with row and column names set accordingly.
# A filepath to a tab-delimited file. This file should contain a gene expression matrix with genes as rows and cells as columns. The first row must contain the cell IDs (header), and the first column must contain the gene names that can have a column name or an empty header.
# A Seurat object. This object should contain gene expression values to be used, stored in object[["RNA"]]@slot_type, where slot_type is the name of the assay slot containing the gene expression matrix to use for prediction (can be either counts or data).
# A filepath to an .rds file containing a Seurat object. This file should contain a Seurat object as described in option 3.
# Please make sure that the input does not contain duplicate gene names or cell IDs.
cytotrace2_res <- cytotrace2(scobj, #seurat对象
                             is_seurat = TRUE, 
                             slot_type = "counts", #counts和data都可以
                             species = 'human')#物种要选择，默认是小鼠
class(cytotrace2_res)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

annotation <- scobj@meta.data

# plotting-一次性生成多个图，然后储存在一个list，用$查看即可
plots <- plotData(cytotrace2_result = cytotrace2_res, 
                  annotation = annotation %>% dplyr::select(cell_type), 
                  is_seurat = TRUE)

plots$CytoTRACE2_UMAP
# ggsave("CytoTRACE2_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Potency_UMAP
# ggsave("CytoTRACE2_Potency_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Relative_UMAP
# ggsave("CytoTRACE2_Relative_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$Phenotype_UMAP
# ggsave("Phenotype_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Boxplot_byPheno
# ggsave("CytoTRACE2_Boxplot_byPheno.pdf",  width = 9, height = 7, dpi = 300)

save(cytotrace2_res, plots, file = "scobj cytotrace2.Rdata")
