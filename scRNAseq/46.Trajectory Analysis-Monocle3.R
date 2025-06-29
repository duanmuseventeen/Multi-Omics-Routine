require(Seurat)
require(monocle3)

cds <- as.cell_data_set(seurat.obj, assay = "RNA")
# 警告: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object

cds <- cluster_cells(cds)
# ·k   Integer number of nearest neighbors to use when creating the k nearest neighbor graph for Louvain/Leiden clustering. k is related to the resolution of the clustering result, a bigger k will result in lower resolution and vice versa. Default is 20.
# ·cluster_method	String indicating the clustering method to use. Options are "louvain" or "leiden". Default is "leiden". Resolution parameter is ignored if set to "louvain".

# # https://github.com/satijalab/seurat-wrappers/issues/114
# integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
# cds <- as.cell_data_set(integrated.sub)

# setting close_loop = FALSE to make result reasonable and simple to understand
# https://github.com/cole-trapnell-lab/monocle3/issues/130
# use_partition = FALSE when the number of partition > 1
cds <- learn_graph(cds, close_loop = FALSE) 
plot_cells(cds,
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)
plot_cells(cds,
           color_cells_by = "cell_type",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

get_earliest_principal_node <- function(cds, cell_type = "Tn"){
  cell_ids <- which(colData(cds)[, "cell_type"] == cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[
      as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
    ]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

all(colnames(seurat.obj@assays$RNA$counts) == names(pseudotime(cds)))
# [1] TRUE
seurat.obj[['pseudotime']] <- pseudotime(cds)
seurat.obj@reductions$umap@cell.embeddings %>% 
  as.data.frame %>% 
  tibble::rownames_to_column("barcode") %>% 
  left_join(seurat.obj@meta.data %>% 
              tibble::rownames_to_column("barcode"), by = "barcode") %>% 
ggplot(aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point() +
  scale_color_viridis(option = "A")

# Working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3, reduction_method = "UMAP")
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# https://stackoverflow.com/questions/44048347/r-open-plotly-in-standalone-window
# cds <- new_cell_data_set(seurat.obj@assays$RNA$counts,
#                          cell_metadata = seurat.obj@meta.data,
#                          gene_metadata = data.frame(
#                            gene_short_name = seurat.obj@assays$RNA$counts,
#                            row.names = seurat.obj@assays$RNA$counts,
#                            stringsAsFactors = FALSE
#                          ))

# https://github.com/rstudio/rstudio/issues/14603
print_app <- function(widget, output_path = "G:/seurat.obj monocle3 3D") {
  
  # Generate random file name
  temp <- paste(output_path, 'html', sep = '.')
  
  # Save. Note, leaving selfcontained=TRUE created files that froze my browser
  htmlwidgets::saveWidget(widget, temp, selfcontained = FALSE)
  
  # Launch with desired application
  system(sprintf("chromium-browser -app=file://%s", temp))
  
  # Return file name if it's needed for any other purpose
  temp
}
print_app(cds_3d_plot_obj)

save(cds, cds_3d, file = "monocle3.Rdata")

rowData(cds) <- data.frame(
  id = rownames(rowData(cds)),
  gene_short_name = rownames(rowData(cds)),
  num_cells_expressed = rowSums(exprs(cds) > 0),
  row.names = rownames(rowData(cds)),
  stringsAsFactors = FALSE
)
cds <- estimate_size_factors(cds)
cds <- estimate_size_factors(cds)
cds <- monocle3::detect_genes(cds)
geneexpr <- fData(cds)
cds <- cds[rownames(rowData(cds)) %in% 
            names(geneexpr$num_cells_expressed)[geneexpr$num_cells_expressed > 10],]

# https://github.com/cole-trapnell-lab/monocle3/issues/344
summary(pseudotime(cds))
cds <- cds[, pseudotime(cds) != Inf]
summary(pseudotime(cds))

gene_fits1 <- fit_models(cds, model_formula_str = "~pseudotime") 
gene_fits2 <- fit_models(cds, model_formula_str = "~pseudotime + orig.ident")
compare_models(gene_fits1, gene_fits2) %>% select(gene_short_name, q_value)
# Likelihood based analysis and quasipoisson
# The quasi-poisson distribution doesn't have a real likelihood function, so some of Monocle's methods won't work with it. Several of the columns in results tables from evaluate_fits() and compare_models() will be NA.

fit_coefs <- coefficient_table(gene_fits2)
fit_coefs <- fit_coefs %>% filter(term != "(Intercept)")
fit_coefs %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
plot_genes_violin(cds, group_cells_by="cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
# plot_genes_hybrid(cds, group_cells_by="cell_type", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))

# expression_family	Distribution	Accuracy	Speed	Notes
# quasipoisson	Quasi-poisson	++	++	Default for fit_models(). Recommended for most users.
# negbinomial	Negative binomial	+++	+	Recommended for users with small datasets (fewer than 1,000 cells).
# poisson	Poisson	-	+++	Not recommended. For debugging and testing only.
# binomial	Binomial	++	++	Recommended for single-cell ATAC-seq

# https://github.com/satijalab/seurat-wrappers/issues/97
AFD_genes = c("SELL","IL7R","HAVCR2","PDCD1","IFNG","GNLY","PRF1","CD44")
cds = cds[rownames(rowData(cds)) %in% AFD_genes,]
cds <- detect_genes(cds)
cds <- cds[, pseudotime(cds) != Inf]
monocle3::plot_genes_in_pseudotime(
  cds,
  label_by_short_name = FALSE,
  color_cells_by="cell_type",
  nrow = 4, ncol = 2,
  min_expr=0.5)
# nrow
# ncol
# panel_order
# vertical_jitter
# horizontal_jitter
# key:   
# model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
# model_expectation <- model_predictions(model_tbl, new_data = new_data)
if(FALSE){
  all(colnames(seurat.obj@assays$RNA$counts) == names(pseudotime(cds)))
  # [1] TRUE
  seurat.obj[['pseudotime']] <- pseudotime(cds)
  dat.gg <- seurat.obj@assays$RNA@counts[AFD_genes,] %>% t %>% 
    as.data.frame %>%
    tibble::rownames_to_column("barcode") %>% 
    left_join(seurat.obj@meta.data %>% 
                tibble::rownames_to_column("barcode"), by = "barcode") %>% 
    tidyr::pivot_longer(cols = AFD_genes, names_to = "Symbol", values_to = "counts")
  ggplot(dat.gg %>% filter(Symbol %in% c("CXCL13")), 
         aes(x = pseudotime, y = counts, color = pseudotime)) +
    geom_point() +
    geom_smooth(method = "loess", n = 50, span = 0.2) +
    scale_color_viridis(option = "A") +
    facet_grid(Symbol ~ .) +
    theme_bw() +
    theme(text = element_text(size = 12), panel.grid = element_blank())
}
