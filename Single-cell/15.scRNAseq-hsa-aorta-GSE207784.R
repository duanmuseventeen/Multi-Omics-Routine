# Aortic Cellular Diversity and Quantitative GWAS Trait Prioritization through Single Nuclear RNA Sequencing (snRNA-Seq) of the Aneurysmal Human Aorta
# GSE207784

# Single nucleus RNA sequencing data are publicly available at the Broad 
# Institute’s Single Cell Portal (https://singlecell.broadinstitute.org/single_cell). 
# The raw dataset is available at the National Center for Biotechnology Information’s 
# Gene Expression Omnibus Database (accession #GSE207784). All other data are contained 
# within the article and its supplementary information, or are available upon reasonable 
# request to the corresponding author.Single nucleus RNA sequencing data are publicly 
# available at the Broad Institute’s Single Cell Portal 
# (https://singlecell.broadinstitute.org/single_cell). The raw dataset is available 
# at the National Center for Biotechnology Information’s Gene Expression Omnibus 
# Database (accession #GSE207784). All other data are contained within the article 
# and its supplementary information, or are available upon reasonable request to 
# the corresponding author.

# We obtained ascobjnding aortic tissue from 6 patients with and 7 patients without 
# thoracic aortic aneurysm and normal tricuspid valves. Baseline characteristics of 
# participants are available in Supplemental Table 1. 

rm(list = ls());gc()
# R----
if(TRUE){
  # Load pkgs---------------------------------------------------------------------
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(clustree) # use for determine the optimal resolution
  library(ROGUE) # use for determine the optimal resolution
  library(harmony)
  # Load Data---------------------------------------------------------------------
  scobj <- "scRNAseq/" %>% 
    Read10X %>% 
    CreateSeuratObject(
      project = "GSE207784",
      min.cells = 0,
      min.features = 0)
  # Add meta.data-----------------------------------------------------------------
  meta <- data.table::fread("AorticAneurysm_MetaData_V1.txt") %>%
    as.data.frame 
  meta.data <- scobj@meta.data %>%
    tibble::rownames_to_column("NAME") %>% 
    left_join(meta, by = "NAME")
  rownames(meta.data) <- meta.data$NAME
  
  scobj@meta.data <- meta.data
  # QC----------------------------------------------------------------------------
  scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, assay = "RNA", pattern = "^MT-")
  scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, assay = "RNA", pattern = "^RP[SL]")
  
  pdf("P2-QC.pdf")
  VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 4)
  FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
  data.frame(
    nFeature = scobj@meta.data$nFeature_RNA
  ) %>% 
    ggplot(aes(x = nFeature)) +
    geom_density() + 
    scale_x_continuous(breaks = c(200,500,1000,3000,4000,5000))
  dev.off()
  
  # In addition to removing clusters of low-quality nuclei with increased cytoplasmic 
  # material, we also performed a per-cluster, per-sample quality control step. In brief, 
  # for each cluster/sample combination, we removed nuclei with extremely low, first 
  # quartile (Q1) minus 1.5 times the interquartile range (IQR), or extremely high, 
  # third quartile (Q3) plus 1.5 times IQR, values for the following metrics: 
  #   1. number of UMI (nUMI), 2. number of unique genes (nGenes), 3. entropy, and 
  # 4. log(nGenes)*entropy. Similarly, we removed nuclei with extremely high values, 
  # Q3 + 1.5 times IQR, for: 1. exon_prop, 2. percent_mito, and 3. doublet_score. 
  # Additional hard cutoffs for removal were set at nGenes < 150, nUMI < 150, or 
  # percent_mito > 5%. When a given cluster/sample combination had less than 30 nuclei, 
  # we removed nuclei with nUMI > 15000 or nUMI ≤ 150, nGenes > 6000 or nGenes ≤ 150, 
  # percent_mito > 5%, entropy < 8, exon_prop > 0.18, doublet_score > 0.3, or 
  # log(nGenes)*entropy > 75.
  
  scobj <- subset(
    scobj, 
    subset = nFeature_RNA > 150 & 
      nFeature_RNA < 6000 & 
      nCount_RNA < 15000 &
      nCount_RNA >= 150 &
      # entropy >= 8 &
      # exon_prop <= 0.18 &
      # doublet_score <= 0.3 &
      # log(nGenes)*entropy > 75 &
      percent.mt < 5)  
  
  pdf("P2-QC-2.pdf")
  VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 4)
  dev.off()
  
  save(scobj, file = "GSE207784-scobj.Rdata")
  
  # After subsetting our data to the top 2000 most highly variable genes, we scaled the data to unit variance and zero mean for each gene and performed principal component (PC) analysis to estimate the first 50 PCs.
  scobj <- scobj %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50)
  
  #Examine and visualize PCA results a few different ways
  pdf("QC-PCA.pdf")
  DimHeatmap(scobj, dims = 1:6, cells = 500, balanced = TRUE)
  ElbowPlot(scobj, reduction = "pca", ndims = 50)
  dev.off()
  # Without Integration-----------------------------------------------------------
  if(TRUE){
    # Principal-component (PC) analysis was performed   on the 2,000 most variable genes, and the first 20 PCs were used for t-SNE and UMAP for data embedding into two dimensions.
    
    # Cluster---------------------------------------------------------------------
    pdf("clustree without integrating data.pdf")
    scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30)
    scobj <- FindClusters(scobj, resolution = seq(0.1, 1, 0.1))
    clustree(scobj, prefix = "RNA_snn_res.")
    dev.off()
    
    # Visualization---------------------------------------------------------------
    #UMAP
    scobj <- RunUMAP(scobj, reduction = "pca", min_dist = 0.3, dims = 1:20)
    #T-SNE
    scobj <- RunTSNE(scobj, reduction = "pca", dims = 1:20)
    
    pdf("umap and tsne without integrating data.pdf")
    scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.5")
    DimPlot(scobj, group.by = "donor_id", reduction = "umap", label = T)
    DimPlot(scobj, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)
    DimPlot(scobj, group.by = "donor_id", reduction = "tsne", label = T)
    DimPlot(scobj, group.by = "RNA_snn_res.0.5", reduction = "tsne", label = T)
    dev.off()
  }
  # Integration with harmony-------------------------------------------------------------------
  if (TRUE) {
    # Cluster---------------------------------------------------------------------
    # To account for biological heterogeneity between samples and technical differences between libraries, we adjusted these PCs using Harmony as implemented in harmony-pytorch v0.1.4 (https://github.com/lilab-bcb/harmony-pytorch) with each sample as their own batch.28  
    scobj.harmony <- RunHarmony(
      scobj,
      group.by.vars = "donor_id",
      reduction.use = "pca",
      reduction.save = "harmony")
    
    # We then constructed a neighborhood graph on these adjusted PCs with 15 neighbors using Euclidean distance followed by Uniform Manifold Approximation and Projection (UMAP) for visualization with min_dist=0.2. Leiden clustering at high resolution (2.0) was used to over-cluster the nuclei into 51 groups.
    scobj.harmony <- FindNeighbors(scobj.harmony, reduction = "harmony", k.param = 15)
    scobj.harmony <- FindClusters(scobj.harmony, resolution = c(seq(0.1, 1, 0.1), 2))
    
    pdf("clustree with integrating data.pdf")
    clustree(scobj.harmony, prefix = "RNA_snn_res.")
    dev.off()
    
    # Visualization---------------------------------------------------------------
    #UMAP
    scobj.harmony <- RunUMAP(scobj.harmony, reduction = "harmony", min_dist = 0.2, dims = 1:10)
    #T-SNE
    scobj.harmony <- RunTSNE(scobj.harmony, reduction = "harmony", dims = 1:10)
    
    pdf("umap and tsne with integrating data.pdf")
    scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.1")
    DimPlot(scobj.harmony, group.by = "donor_id", reduction = "umap", label = T)
    DimPlot(scobj.harmony, group.by = "RNA_snn_res.1", reduction = "umap", label = T)
    DimPlot(scobj.harmony, group.by = "donor_id", reduction = "tsne", label = T)
    DimPlot(scobj.harmony, group.by = "RNA_snn_res.1", reduction = "tsne", label = T)
    dev.off()
    
    save(scobj.harmony, file = "GSE207784-scobj.harmony.Rdata")
    # Annotation----------------------------------------------------------------
    load("GSE207784-scobj.harmony.Rdata")
    scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.6")
    
    # The most common cell type observed in aortas were VSMCs, which comprised 
    # approximately ~65% of nuclei in our samples (clusters 1, 2, 3), consistent 
    # with observations from aortic histology. (Figure 2c, Supplemental Table 3) 
    # Clusters 1, 2, 3, and 4 all expressed canonical markers of VSMCs including 
    # MYH11 and ACTA2, and included VSMC subtypes as well as pericytes.
    
    # Fibroblast cluster 6 was the next most common cell type after VSMCs and expressed typical but less specific gene markers such as LAMA2 and TSHZ2.
    
    # As expected for vascular tissue, endothelial cells expressing PECAM1 (CD31) were also numerous and were observed in three separate clusters 8, 9, 10.
    
    SMC <- c('MYH11', 'ACTA2','CNN1', 'SMTN', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # ?0 1 3 ?5 6 ?7 9 ?12 
    Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 3
    FB <- c('LAMA2', 'TSHZ2') # 2 
    EC <- c('PECAM1') # 4
    RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ') # ?12
    Meso <- c('MSLN', 'WT1', 'BNC1') # 7
    Adi <- c('GPAM', 'FASN', 'LEP') # 
    Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 
    Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 
    M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74', 'CD14', 'C1QA') # 
    
    VlnPlot(scobj.harmony, features = M, group.by = "RNA_snn_res.0.6", layer = "counts", log = TRUE)
    
    scobj.harmony@meta.data$cell_type <- "Unknown"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(0,1,5,6,7,9)] <- "SMC"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(3)] <- "Peri"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(7)] <- "Meso"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(2)] <- "FB"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(4)] <- "EC"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(8)] <- "Macrophage"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(10)] <- "Lymphoid"
    scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.6 %in% c(12)] <- "Neural"
    
    p1 <- DimPlot(scobj.harmony, reduction = "umap",
                  group.by = "RNA_snn_res.0.6",
                  label = TRUE, pt.size = 0.5) 
    p2 <- DimPlot(scobj.harmony, reduction = "umap",
                  group.by = "cell_type",
                  label = TRUE, pt.size = 0.5)
    p3 <- DimPlot(scobj.harmony, reduction = "umap",
                  group.by = "donor_id",
                  label = TRUE, pt.size = 0.5) 
    p1 + p2 + p3
    # ggsave(..., width = 18, height = 6)
    
    FeaturePlot(object = scobj.harmony, slot = "counts", features = 'Targetgene')
    RidgePlot(object = scobj.harmony, features = 'Targetgene')
    
    save(scobj.harmony, file = "GSE207784-scobj.annot.Rdata")
  }
  # Integration with cca-------------------------------------------------------------------
  if (TRUE) {
    # Cluster---------------------------------------------------------------------
    # SplitObject(scobj, split.by = "donor_id") # into list
    scobj.split <- split(scobj, f = scobj@meta.data$donor_id %>% factor) # split RNA layer
    
    scobj.split <- scobj.split %>% 
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
      ScaleData %>% 
      RunPCA(npcs = 50)
    
    # To account for biological heterogeneity between samples and technical differences between libraries, we adjusted these PCs using Harmony as implemented in harmony-pytorch v0.1.4 (https://github.com/lilab-bcb/harmony-pytorch) with each sample as their own batch.28  
    scobj.cca <- IntegrateLayers(
      object = scobj.split, 
      method = CCAIntegration,
      orig.reduction = "pca", 
      new.reduction = "integrated.cca",
      verbose = FALSE
    )
    
    # re-join layers after integration
    scobj.cca[["RNA"]] <- JoinLayers(scobj.cca[["RNA"]])
    
    # We then constructed a neighborhood graph on these adjusted PCs with 15 neighbors using Euclidean distance followed by Uniform Manifold Approximation and Projection (UMAP) for visualization with min_dist=0.2. Leiden clustering at high resolution (2.0) was used to over-cluster the nuclei into 51 groups.
    scobj.cca <- FindNeighbors(scobj.cca, reduction = "integrated.cca", k.param = 15)
    scobj.cca <- FindClusters(scobj.cca, resolution = c(seq(0.1,1,0.1), 2))
    
    pdf("clustree with integrating cca.pdf")
    clustree(scobj.cca, prefix = "RNA_snn_res.")
    dev.off()
    
    # Visualization---------------------------------------------------------------
    #UMAP
    scobj.cca <- RunUMAP(scobj.cca, reduction = "integrated.cca", min_dist = 0.2, dims = 1:10)
    #T-SNE
    scobj.cca <- RunTSNE(scobj.cca, reduction = "integrated.cca", dims = 1:10)
    
    pdf("umap and tsne with integrating cca.pdf")
    scobj.cca <- RegroupIdents(scobj.cca, metadata = "RNA_snn_res.1")
    DimPlot(scobj.cca, group.by = "donor_id", reduction = "umap", label = T)
    DimPlot(scobj.cca, group.by = "RNA_snn_res.1", reduction = "umap", label = T)
    DimPlot(scobj.cca, group.by = "donor_id", reduction = "tsne", label = T)
    DimPlot(scobj.cca, group.by = "RNA_snn_res.1", reduction = "tsne", label = T)
    dev.off()
  }
  # Annotation--------------------------------------------------------------------
  # the integration effect of cca is worse than harmony according to clustree
  # Omit
}
# Python----
if(TRUE){
# https://blog.csdn.net/u011262253/article/details/118733569

# Load pkgs---------------------------------------------------------------------
import numpy as np
import pandas as pd
import scanpy as sc
import pooch
from harmony import harmonize

# Settings----------------------------------------------------------------------
sc.settings.set_figure_params(dpi = 150, facecolor='white')
# 用于设置日志的详细程度
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

# Load data---------------------------------------------------------------------
# scobj = sc.read_10x_mtx(
#   'scRNAseq/',  # 不支持压缩格式，文件需提前解压
#   var_names='gene_symbols',             # 使用 gene_symbols 作为变量名
#   cache=True)                           # 写入缓存，可以更快的读取文件
  
# https://www.jianshu.com/p/8c970e2c427a
adata = sc.read_mtx('scRNAseq/matrix.mtx')
adata_bc=pd.read_csv('scRNAseq/barcodes.tsv',header=None)
adata_features=pd.read_csv('scRNAseq/features.tsv',header=None,sep="\t")
adata= adata.T
adata.obs= adata_bc
adata.obs.columns = ['cell_id']
adata.var['gene_name']= adata_features[1].tolist()
adata.var.index= adata.var['gene_name']
# Add meta.data-----------------------------------------------------------------
# https://www.jianshu.com/p/8627435e9414

metadata=pd.read_excel("MetaData.xlsx")

all(adata.obs["cell_id"] == metadata["cell_id"])

adata.obs= metadata
# QC----------------------------------------------------------------------------
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
  
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
  adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
  adata,
  ["n_genes_by_counts", "total_counts", "pct_counts_mt","pct_counts_ribo","pct_counts_hb"],
  jitter=0.4,
  multi_panel=True,
)

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")

# We then log-normalized the count data by dividing the UMI for each gene in a droplet by the total number of UMI in the droplet, multiplying by 10,000, and taking the natural log.
sc.pp.log1p(adata)

# After subsetting our data to the top 2000 most highly variable genes, 
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]

# we scaled the data to unit variance and zero mean for each gene
sc.pp.regress_out(adata, ["total_counts"])
sc.pp.scale(adata)

# performed principal component (PC) analysis to estimate the first 50 PCs
sc.tl.pca(adata, svd_solver="arpack")

sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

adata.obs['donor_id'] = pd.Categorical(adata.obs['donor_id'])
adata.obs['disease'] = pd.Categorical(adata.obs['disease'])
sc.pl.pca(
  adata,
  color=["donor_id", "disease", "pct_counts_mt", "pct_counts_ribo"],
  dimensions=[(0, 1), (0, 1), (0, 1), (0, 1)],
  ncols=2,
  size=2,
)

# we adjusted these PCs using Harmony as implemented in harmony-pytorch v0.1.4 (https://github.com/lilab-bcb/harmony-pytorch) with each sample as their own batch
Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'donor_id')
adata.obsm['X_harmony'] = Z

# We then constructed a neighborhood graph on these adjusted PCs with 15 neighbors using Euclidean distance followed by Uniform Manifold Approximation and Projection (UMAP) for visualization with min_dist=0.2.
sc.pp.neighbors(adata, n_neighbors=15, metric='euclidean')
sc.tl.umap(adata, min_dist = 0.2)

# Leiden clustering at high resolution (2.0) was used to over-cluster the nuclei into 51 groups.
for res in [0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0]:
  sc.tl.leiden(
    adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
  )

sc.pl.umap(
  adata,
  color=["leiden_res_0.20", "leiden_res_0.50", "leiden_res_1.00" , "leiden_res_2.00"],
  legend_loc="on data",
)

sc.pl.umap(
  adata,
  color=["leiden_res_0.50", "leiden_res_1.00" , "biosample_id"],
  legend_loc="on data",
)

# the batch effect is significant after integration
}






