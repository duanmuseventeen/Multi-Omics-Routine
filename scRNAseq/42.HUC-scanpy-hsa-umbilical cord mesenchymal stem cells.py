# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
# https://scanpy.readthedocs.io/en/stable/tutorials/basics/integrating-data-using-ingest.html
# https://github.com/brianhie/scanorama

# Load pkgs--------------------------------------------------------------------
import os
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
# Load data--------------------------------------------------------------------
# sc._settings.ScanpyConfig.autosave

path = 'rawdata/'
print(os.listdir(path))

genes = None
cell = None
mtx = None
scobj_list = {}

for name in os.listdir(path):
    scobj = sc.read_10x_mtx(path + name, var_names='gene_symbols')
    scobj.var_names_make_unique()
    scobj_list[name] = scobj

scobj_list

HUC = ad.concat(scobj_list, label="sample")
# UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
#   utils.warn_names_duplicates("obs")
HUC.obs_names_make_unique()
print(HUC.obs["sample"].value_counts())
HUC

# sample
# hUC_1    8778
# hUC_2    8333
# hUC_3    7718

HUC.shape
# (24829, 36601)
# Quality Control--------------------------------------------------------------------
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
HUC.var["mt"] = HUC.var_names.str.startswith("MT-")
# ribosomal genes
HUC.var["ribo"] = HUC.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
HUC.var["hb"] = HUC.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    HUC, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    HUC,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt",
        "pct_counts_ribo", "pct_counts_hb"],
    jitter=0.4,
    multi_panel=True,
    groupby="sample",
    save=" QC.pdf"
)
sc.pl.scatter(HUC, "total_counts", "n_genes_by_counts",
              color="sample", save=" QC.pdf")


# sc.pp.filter_genes(adata, max_counts=15000)
sc.pp.filter_cells(HUC, min_genes=1000)
sc.pp.filter_cells(HUC, max_genes=5000)
sc.pp.filter_cells(HUC, max_counts=15000)

# 获取线粒体基因占比在 5% 以下的细胞样本
HUC = HUC[HUC.obs.pct_counts_mt < 5, :]
# 获取红细胞基因数在 1% 以下的细胞样本
HUC = HUC[HUC.obs.pct_counts_hb < 1, :]

HUC.shape
# (17683, 36601)

sc.pl.violin(
    HUC,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt",
        "pct_counts_ribo", "pct_counts_hb"],
    jitter=0.4,
    multi_panel=True,
    groupby="sample",
    save=" after QC.pdf"
)
# Doublet detection--------------------------------------------------------------------
# Alternative methods for doublet detection within the scverse ecosystem are DoubletDetection and SOLO.
# You can read more about these in the Doublet Detection chapter of Single Cell Best Practices.
sc.pp.scrublet(HUC, batch_key="sample")

singlet = HUC.obs['predicted_doublet']
singlet.value_counts()
# predicted_doublet
# False    17668
# True        15

HUC = HUC[HUC.obs.predicted_doublet == False, :]
HUC.shape
# (17668, 36601)
# Normalization--------------------------------------------------------------------
# Saving count data
HUC.layers["counts"] = HUC.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(HUC)
# Logarithmize the data
sc.pp.log1p(HUC)
# Feature selection--------------------------------------------------------------------
sc.pp.highly_variable_genes(HUC, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(HUC)
# Dimensionality Reduction--------------------------------------------------------------------
# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.

sc.tl.pca(HUC)

sc.pl.pca_variance_ratio(HUC, n_pcs=50, log=True)

sc.pl.pca(
    HUC,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
# Batch effect(if necessary)--------------------------------------------------------------------
# Scanorama [Hie et al., 2019] is an algorithm for integrating single-cell data from multiple experiments stored in an AnnData object.
# This function should be run after performing PCA but before computing the neighbor graph, as illustrated in the example below.

sce.pp.scanorama_integrate(HUC, 'sample', verbose=1)
# Nearest neighbor graph constuction and visualization--------------------------------------------------------------------
sc.pp.neighbors(HUC)

sc.tl.umap(HUC)

sc.pl.umap(
    HUC,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
    save="2000 umap.pdf"
)

# If you inspect batch effects in your UMAP it can be beneficial to integrate across samples and
# perform batch correction/integration. We recommend checking out scanorama and scvi-tools for
# batch integration.
# Clustering--------------------------------------------------------------------
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(HUC, flavor="igraph", n_iterations=2)
sc.pl.umap(HUC, color=["leiden"], save="cluster.pdf")

# BBKNN-------------------------------------------------------------------------
# sc.external.pp.bbknn(HUC, batch_key="batch")
# sc.tl.umap(HUC)
# Differentially-expressed Genes as Markers--------------------------------------------------------------------
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(HUC, groupby="leiden_res_0.50", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    HUC, groupby="leiden_res_0.50", standard_scale="var", n_genes=5
)

sc.get.rank_genes_groups_df(HUC, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(HUC, group="7").head(5)["names"]
sc.pl.umap(
    HUC,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)

# You may have noticed that the p-values found here are extremely low. This is due to the
# statistical test being performed considering each cell as an independent sample.
# For a more conservative approach you may want to consider “pseudo-bulking” your data by sample
# (e.g. sc.get.aggregate(adata, by=["sample", "cell_type"], func="sum", layer="counts")) and
# using a more powerful differential expression tool, like pydeseq2.
# Manual cell-type annotation--------------------------------------------------------------------
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "COL4A4",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
    ],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}
