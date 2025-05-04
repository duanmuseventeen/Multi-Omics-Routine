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