library(monocle)
library(Seurat)

# Note: Don't load monocle and monocle3 together

#============================== Monocle2========================================
# REF:
# https://www.jianshu.com/p/5d6fd4561bc0
# https://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle

pd <- new("AnnotatedDataFrame", data = macro@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(
  id = rownames(macro),
  gene_short_name = rownames(macro),
  num_cells_expressed = rowSums(macro@assays$RNA$counts > 0),
  row.names = rownames(macro),
  stringsAsFactors = FALSE
))
cds2 <- monocle::newCellDataSet(
  cellData = macro@assays$RNA$counts, 
  phenoData = pd, 
  featureData = fd,
  expressionFamily=negbinomial.size())

cds2 <- estimateSizeFactors(cds2)
cds2 <- estimateDispersions(cds2)
# cds2 <- monocle::clusterCells(cds2)
# 错误于monocle::clusterCells(cds2):
#   k must be smaller than the total number of points!
# plot_cell_clusters(cds2, 1, 2, color = "cell_type")

cds2 -> cds2.raw
# Filtering low-quality cells Recommended
cds2 <- detectGenes(cds2, min_expr = 0.1)
print(head(fData(cds2)))
# https://github.com/cole-trapnell-lab/garnett/issues/20
# 错误: 函数‘fData’标签‘x = "CellDataSet"’找不到继承方法
# I suspect that you're right that this is a conflict between monocle2 and monocle3. I have just updated the Garnett website to have separate instructions for those using Monocle2 versus monocle3. Try using one or the other of the instructions and making sure that you don't have both monocle and monocle3 loaded.
expressed_genes <- row.names(subset(fData(cds2), num_cells_expressed >= 10))

pData(cds2)$Total_mRNAs <- Matrix::colSums(exprs(cds))
# cds <- cds[,pData(cds)$Total_mRNAs < 1e6]
# all(pData(cds)$Total_mRNAs < 1e6)
# [1] TRUE

upper_bound <- 10^(mean(log10(pData(cds2)$Total_mRNAs)) +
                     2*sd(log10(pData(cds2)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds2)$Total_mRNAs)) -
                     2*sd(log10(pData(cds2)$Total_mRNAs)))

# qplot(Total_mRNAs, data = pData(cds), color = cell_type, geom =
#         "density") +
#   geom_vline(xintercept = lower_bound) +
#   geom_vline(xintercept = upper_bound)

cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
             pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)

# # Log-transform each value in the expression matrix.
# L <- log(exprs(HSMM[expressed_genes,]))
# 
# # Standardize each gene, so that they are all on the same scale,
# # Then melt the data with plyr so we can plot it easily
# melted_dens_df <- reshape2::melt(Matrix::t(scale(Matrix::t(L))))
# 
# # Plot the distribution of the standardized gene expression values.
# qplot(value, geom = "density", data = melted_dens_df) +
#   stat_function(fun = dnorm, size = 0.5, color = 'red') +
#   xlab("Standardized log(FPKM)") +
#   ylab("Density")

# Constructing Trajectory
diff_test_res <- differentialGeneTest(
  cds[expressed_genes,],
  fullModelFormulaStr = "~sm.ns(Pseudotime)")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds, root_state = "Tn_CCR7")
plot_cell_trajectory(cds, color_by = "cell_type")
plot_cell_trajectory(cds, color_by = "State")

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))
plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

blast_genes <- row.names(subset(fData(cds), 
                                gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)

cds_expressed_genes <-  row.names(subset(fData(cds),
                                         num_cells_expressed >= 10))
cds_filtered <- cds[cds_expressed_genes,]
my_genes <- row.names(subset(fData(cds_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- cds_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
