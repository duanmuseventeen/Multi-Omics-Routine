# Reference
# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html

# Preprocessing and clustering 3k PBMCs (legacy workflow)-----------------------
# In May 2017, this started out as a demonstration that Scanpy would allow to reproduce most of Seurat’s guided clustering tutorial (Satija et al., 2015).

# We gratefully acknowledge Seurat’s authors for the tutorial! In the meanwhile, we have added and removed a few pieces.

# The data consist of 3k PBMCs from a Healthy Donor and are freely available from 10x Genomics (here from this webpage). On a unix system, you can uncomment and run the following to download and unpack the data. The last line creates a directory for writing processed data.

# Linux-------------------------------------------------------------------------
# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write
# Load pkgs---------------------------------------------------------------------
import pandas as pd
import scanpy as sc

# 用于设置日志的详细程度
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

results_file = "write/pbmc3k.h5ad"  # the file that will store the analysis results

# Read in the count matrix into an AnnData object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5-based file format: .h5ad.
# 使用绝对路径，相对路径报错
adata = sc.read_10x_mtx(
  "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
  var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
  cache=True,  # write a cache file for faster subsequent reading
)

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
adata
# AnnData object with n_obs × n_vars = 2700 × 32738
# var: 'gene_ids'
# Preprocessing-----------------------------------------------------------------
# Show those genes that yield the highest fraction of counts in each single cell, across all cells.

sc.pl.highest_expr_genes(adata, n_top=20)
# normalizing counts per cell
# finished (0:00:00)

# Basic filtering:
  
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# filtered out 19024 genes that are detected in less than 3 cells
# Let’s assemble some information about mitochondrial genes, which are important for quality control.

# Citing from “Simple Single Cell” workflows (Lun, McCarthy & Marioni, 2017):
#   
# High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.
# 
# With pp.calculate_qc_metrics, we can compute many metrics very efficiently.

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
  adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
# A violin plot of some of the computed quality measures:
#   
# · the number of genes expressed in the count matrix
# · the total counts per cell
# · the percentage of counts in mitochondrial genes

sc.pl.violin(
  adata,
  ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
  jitter=0.4,
  multi_panel=True,
)
# Remove cells that have too many mitochondrial genes expressed or too many total counts:
  
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
# Actually do the filtering by slicing the AnnData object.

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
# Total-count normalize (library-size correct) the data matrix 
# to 10,000 reads per cell, so that counts become comparable among cells.

sc.pp.normalize_total(adata, target_sum=1e4)
# normalizing counts per cell
# finished (0:00:00)
# Logarithmize the data:
sc.pp.log1p(adata)
  
# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)
# Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

# Note
# You can get back an AnnData of the object in .raw by calling .raw.to_adata().

adata.raw = adata

# Note
# If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
# The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.

# Actually do the filtering

adata = adata[:, adata.var.highly_variable]
Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)
# Principal component analysis--------------------------------------------------
# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.

sc.tl.pca(adata, svd_solver="arpack")

# We can make a scatter plot in the PCA coordinates, but we will not use that later on.

sc.pl.pca(adata, color="CST3")

# Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function sc.tl.louvain() or tSNE sc.tl.tsne(). In our experience, often a rough estimate of the number of PCs does fine.
sc.pl.pca_variance_ratio(adata, log=True)

adata.write(results_file)
adata
# AnnData object with n_obs × n_vars = 2638 × 1838
# obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
# var: 'gene_ids', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'
# uns: 'log1p', 'hvg', 'pca'
# obsm: 'X_pca'
# varm: 'PCs'

# Computing the neighborhood graph----------------------------------------------
# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# computing neighbors
# using 'X_pca' with n_pcs = 40
# finished: added to `.uns['neighbors']`
# `.obsp['distances']`, distances for each pair of neighbors
# `.obsp['connectivities']`, weighted adjacency matrix (0:00:01)

# Embedding the neighborhood graph----------------------------------------------
# We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below. It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
  
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')
sc.tl.umap(adata)
# computing UMAP
# finished: added
# 'X_umap', UMAP coordinates (adata.obsm) (0:00:03)
sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"])

# As we set the .raw attribute of adata, the previous plots showed the “raw” (normalized, logarithmized, but uncorrected) gene expression. You can also plot the scaled and corrected gene expression by explicitly stating that you don’t want to use .raw.

sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"], use_raw=False)
# Clustering the neighborhood graph---------------------------------------------
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag et al. (2018). Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

sc.tl.leiden(
  adata,
  resolution=0.9,
  random_state=0,
  flavor="igraph",
  n_iterations=2,
  directed=False,
)
# running Leiden clustering
# finished: found 8 clusters and added
# 'leiden', the cluster labels (adata.obs, categorical) (0:00:00)

# Plot the clusters, which agree quite well with the result of Seurat.
sc.pl.umap(adata, color=["leiden", "CST3", "NKG7"])

# Save the result.
adata.write(results_file)
# Finding marker genes----------------------------------------------------------
# Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.
sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# ranking genes
# finished: added to `.uns['rank_genes_groups']`
# 'names', sorted np.recarray to be indexed by group ids
# 'scores', sorted np.recarray to be indexed by group ids
# 'logfoldchanges', sorted np.recarray to be indexed by group ids
# 'pvals', sorted np.recarray to be indexed by group ids
# 'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)

sc.settings.verbosity = 2  # reduce the verbosity

# The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar. We recommend using the latter in publications, see e.g., Sonison & Robinson (2018). You might also consider much more powerful differential testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Save the result.
adata.write(results_file)
# As an alternative, let us rank genes using logistic regression. For instance, this has been suggested by Natranos et al. (2018). The essential difference is that here, we use a multi-variate appraoch whereas conventional differential tests are uni-variate. Clark et al. (2014) has more details.

sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# With the exceptions of IL7R, which is only found by the t-test and FCER1A, which is only found by the other two appraoches, all marker genes are recovered in all approaches.

# Let us also define a list of marker genes for later reference.

marker_genes = [
  *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
  *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
  *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]
# Reload the object that has been save with the Wilcoxon Rank-Sum test result.

  adata = sc.read(results_file)
# Show the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.

pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)
# 0	1	2	3	4	5	6	7
# 0	RPS12	CD74	LST1	NKG7	CCL5	LYZ	HLA-DPA1	PF4
# 1	LDHB	CD79A	FCER1G	GZMB	NKG7	S100A9	HLA-DPB1	SDPR
# 2	RPS25	HLA-DRA	AIF1	GNLY	CST7	S100A8	HLA-DRA	GNG11
# 3	RPS27	CD79B	COTL1	CTSW	B2M	TYROBP	HLA-DRB1	PPBP
# 4	RPS6	HLA-DPB1	FCGR3A	PRF1	GZMA	FTL	CD74	NRGN

# Get a table with the scores and groups.

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
pd.DataFrame(
  {
    group + "_" + key[:1]: result[key][group]
    for group in groups
    for key in ["names", "pvals"]
  }
).head(5)
0_n	0_p	1_n	1_p	2_n	2_p	3_n	3_p	4_n	4_p	5_n	5_p	6_n	6_p	7_n	7_p
0	RPS12	4.164150e-226	CD74	2.487145e-183	LST1	4.980059e-88	NKG7	3.591782e-93	CCL5	4.054732e-126	LYZ	2.844372e-249	HLA-DPA1	5.422417e-21	PF4	4.722886e-10
1	LDHB	1.793331e-223	CD79A	1.679730e-170	FCER1G	1.449472e-84	GZMB	2.033412e-87	NKG7	1.528118e-110	S100A9	2.654880e-246	HLA-DPB1	7.591860e-21	SDPR	4.733899e-10
2	RPS25	4.641187e-204	HLA-DRA	6.949695e-167	AIF1	5.839647e-83	GNLY	8.130223e-85	CST7	1.332168e-85	S100A8	8.731315e-238	HLA-DRA	1.306768e-19	GNG11	4.733899e-10
3	RPS27	9.438482e-192	CD79B	2.569135e-154	COTL1	1.261406e-81	CTSW	6.944632e-84	B2M	8.615108e-79	TYROBP	9.799314e-221	HLA-DRB1	1.865104e-19	PPBP	4.744938e-10
4	RPS6	5.780179e-188	HLA-DPB1	3.580735e-148	FCGR3A	4.610698e-77	PRF1	1.621421e-83	GZMA	1.480430e-78	FTL	3.676035e-215	CD74	5.853161e-19	NRGN	4.800511e-10

# Compare to a single cluster:
sc.tl.rank_genes_groups(adata, "leiden", groups=["0"], reference="1", method="wilcoxon")
sc.pl.rank_genes_groups(adata, groups=["0"], n_genes=20)

# If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)

# Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):
adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)

# If you want to compare a certain gene across groups, use the following.
sc.pl.violin(adata, ["CST3", "NKG7", "PPBP"], groupby="leiden")

# Actually mark the cell types.
new_cluster_names = [
  "CD4 T",
  "B",
  "FCGR3A+ Monocytes",
  "NK",
  "CD8 T",
  "CD14+ Monocytes",
  "Dendritic",
  "Megakaryocytes",
]
adata.rename_categories("leiden", new_cluster_names)
sc.pl.umap(
  adata, color="leiden", legend_loc="on data", title="", frameon=False, save=".pdf"
)
# WARNING: saving figure to file figures/umap.pdf

# Now that we annotated the cell types, let us visualize the marker genes.
sc.pl.dotplot(adata, marker_genes, groupby="leiden")
# There is also a very compact violin plot.

sc.pl.stacked_violin(adata, marker_genes, groupby="leiden");
# During the course of this analysis, the AnnData accumlated the following annotations.

adata
# AnnData object with n_obs × n_vars = 2638 × 1838
# obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'leiden'
# var: 'gene_ids', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'
# uns: 'hvg', 'leiden', 'leiden_colors', 'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'umap'
# obsm: 'X_pca', 'X_umap'
# varm: 'PCs'
# obsp: 'connectivities', 'distances'
# `compression='gzip'` saves disk space, and slightly slows down writing and subsequent reading
adata.write(results_file, compression="gzip")
# Get a rough overview of the file using h5ls, which has many options - for more details see here. The file format might still be subject to further optimization in the future. All reading functions will remain backwards-compatible, though.

# If you want to share this file with people who merely want to use it for visualization, a simple way to reduce the file size is by removing the dense scaled and corrected data matrix. The file still contains the raw data used in the visualizations in adata.raw.

adata.raw.to_adata().write("./write/pbmc3k_withoutX.h5ad")
