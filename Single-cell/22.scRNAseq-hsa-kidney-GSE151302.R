# 	Single Cell Transcriptional and Chromatin Accessibility Profiling on the human adult kidneys

# GSE151302

# 一些可能有用的经验：
# 1. subset函数downsample有内置的seed，结果可重复
# 2. 以本数据为例，downsample处理与否，对结果影响不大。downsample后，注释更简单，
#    可以认为downsample处理优于不做downsample
# 3. scdblfinder函数可能需要优化参数设置，如果时间允许，doubletfinder的结果可能更可信

# n=5 human adult kidney cortex samples were analyzed with snRNA-seq and snATAC-seq 
# (10xGenomics platform). snRNA-seq sequencing of 3 samples (Healthy1/2/3) were 
# already done and deposited on GSE131882 as control1/2/3

GSM4572187	Healthy1_snATACseq
GSM4572188	Healthy2_snATACseq
GSM4572189	Healthy3_snATACseq
GSM4572190	Healthy4_snATACseq
GSM4572191	Healthy5_snATACseq
GSM4572192	Healthy1_snRNAseq
GSM4572193	Healthy2_snRNAseq
GSM4572194	Healthy3_snRNAseq
GSM4572195	Healthy4_snRNAseq
GSM4572196	Healthy5_snRNAseq

# Reference---------------------------------------------------------------------
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://satijalab.org/seurat/articles/sctransform_vignette

rm(list = ls());gc()
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(Matrix)
library(scDblFinder)
library(Signac)
library(cicero)
library(sctransform)
library(hdf5r)
library(GenomicRanges)
## snRNAseq---------------------------------------------------------------------
if(TRUE){
  # Load data---------------------------------------------------------------------
  scobjList <- lapply(dir(), function(x){
    tmp <- dir(x)
    sce <- paste0(x, "/", tmp) %>% 
      Read10X_h5 %>% 
      CreateSeuratObject(project = tmp)
    return(sce)
  })
  
  scobj <- merge(x = scobjList[[1]], y = scobjList[-1], add.cell.ids = dir())
  
  head(scobj@meta.data)
  # orig.ident nCount_RNA
  # Control1_AAACCTGAGGGTCTCC-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5       3410
  # Control1_AAACCTGAGTGTTAGA-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5       5341
  # Control1_AAACCTGCAAGCGCTC-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5      10178
  # Control1_AAACCTGCACCAGATT-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5       8911
  # Control1_AAACCTGCACGGCGTT-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5        894
  # Control1_AAACCTGCAGTCAGAG-1 GSM4572192_Control1_filtered_feature_bc_matrix.h5       2035
  # nFeature_RNA
  # Control1_AAACCTGAGGGTCTCC-1         1862
  # Control1_AAACCTGAGTGTTAGA-1         2501
  # Control1_AAACCTGCAAGCGCTC-1         3580
  # Control1_AAACCTGCACCAGATT-1         3018
  # Control1_AAACCTGCACGGCGTT-1          731
  # Control1_AAACCTGCAGTCAGAG-1         1314
  
  table (scobj@meta.data$orig.ident)
  # GSM4572192_Control1_filtered_feature_bc_matrix.h5 
  # 6905 
  # GSM4572193_Control2_filtered_feature_bc_matrix.h5 
  # 4245 
  # GSM4572194_Control3_filtered_feature_bc_matrix.h5 
  # 6599 
  # GSM4572195_Control4_filtered_feature_bc_matrix.h5 
  # 4496 
  # GSM4572196_Control5_filtered_feature_bc_matrix.h5 
  # 4708 
  # QC----------------------------------------------------------------------------
  scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
  scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
  scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")
  
  VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
  FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  data.frame(
    nFeature = scobj@meta.data$nFeature_RNA,
    group = scobj@meta.data$orig.ident
  ) %>% 
    ggplot(aes(x = nFeature, color = group)) +
    geom_density() + 
    scale_x_continuous(breaks = c(200,300,400,500,1000,3000,4000,5000))
  
  # Subsequently datasets were preprocessed with Seurat v3.0.221 to remove low-quality 
  # nuclei (Features > 500, Features < 4000, RNA count < 16000, %Mitochondrial genes < 0.8,
  # %Ribosomal protein large or small subunits < 0.4) and 
  
  scobj <- subset(
    scobj, 
    subset = 
      nFeature_RNA > 500 & 
      nFeature_RNA < 4000 & 
      nCount_RNA < 16000 &
      percent.mt < 5 &
      percent.rp < 2 &
      percent.hb < 1)  
  
  VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
  
  table(scobj@meta.data$orig.ident)
  # GSM4572192_Control1_filtered_feature_bc_matrix.h5 
  # 5689 
  # GSM4572193_Control2_filtered_feature_bc_matrix.h5 
  # 3506 
  # GSM4572194_Control3_filtered_feature_bc_matrix.h5 
  # 5910 
  # GSM4572195_Control4_filtered_feature_bc_matrix.h5 
  # 3948 
  # GSM4572196_Control5_filtered_feature_bc_matrix.h5 
  # 4060 
  
  # scobj.ds <- subset(scobj, downsample = 3500) # by default seed = 1, thus the result is reproducible
  # 
  # table(scobj.ds@meta.data$orig.ident)
  # # GSM4572192_Control1_filtered_feature_bc_matrix.h5 
  # # 3500
  # # GSM4572193_Control2_filtered_feature_bc_matrix.h5 
  # # 3500 
  # # GSM4572194_Control3_filtered_feature_bc_matrix.h5 
  # # 3500 
  # # GSM4572195_Control4_filtered_feature_bc_matrix.h5 
  # # 3500 
  # # GSM4572196_Control5_filtered_feature_bc_matrix.h5 
  # # 3500 
  
  if(TRUE){ 
    # Apply sctransform normalization-----------------------------------------------
    # Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
    # Transformed data will be available in the SCT assay, which is set as the default after running sctransform
    # During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
    # In Seurat v5, SCT v2 is applied by default. You can revert to v1 by setting vst.flavor = 'v1'
    # The glmGamPoi package substantially improves speed and is used by default if installed, with instructions here
    
    # The filtered library was normalized with SCTransform
    
    # While this method of normalization is standard and widely used in scRNA-seq analysis, 
    # global-scaling relies on an assumption that each cell originally contains the same 
    # number of RNA molecules. We and others have developed alternative workflows for the 
    # single cell preprocessing that do not make these assumptions. For users who are 
    # interested, please check out our SCTransform() normalization workflow. The method 
    # is described in ourpaper, with a separate vignette using Seurat here. The use of 
    # SCTransform replaces the need to run NormalizeData, FindVariableFeatures, or 
    # ScaleData (described below.)
    
    scobj.ds.harmony <- scobj %>% 
      SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
      RunPCA(npcs = 50) %>% 
    # Integration-------------------------------------------------------------------
    # corrected for batch effects with Harmony v1.067 using the “RunHarmony” function in Seurat.
    # Perform integration (harmony)
      RunHarmony(
        group.by.vars = "orig.ident",
        reduction.use = "pca",
        reduction.save = "harmony") %>% 
    # Cluster-----------------------------------------------------------------------
      FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
      FindClusters(resolution = seq(0.1, 1, 0.1))
    
    clustree(scobj.ds.harmony, prefix = "SCT_snn_res.")
    # Visualization-----------------------------------------------------------------
    #UMAP
    scobj.ds.harmony <- RunUMAP(scobj.ds.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
    #T-SNE
    scobj.ds.harmony <- RunTSNE(scobj.ds.harmony, reduction = "harmony", dims = 1:30)
    
    p1 <- DimPlot(scobj.ds.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend()
    p2 <- DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.4", reduction = "umap", label = T)
    p1 + p2
    # Cell cycle--------------------------------------------------------------------
    exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
                          header = TRUE, as.is = TRUE, row.names = 1)
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    
    scobj.ds.harmony <- CellCycleScoring(scobj.ds.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    RidgePlot(scobj.ds.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
    
    scobj.ds.harmony <- RunPCA(scobj.ds.harmony, features = c(s.genes, g2m.genes))
    DimPlot(scobj.ds.harmony, reduction = "pca")
    DimPlot(scobj.ds.harmony, reduction = "umap")
    # ScdblFinder---------------------------------------------------------------
    if(TRUE){
      scobj.ds.harmony.sce <- as.SingleCellExperiment(scobj.ds.harmony, assay = "SCT")
      scobj.ds.harmony.sce <- scDblFinder(
        scobj.ds.harmony.sce, 
        samples = scobj.ds.harmony.sce$orig.ident,
        clusters = "SCT_snn_res.0.4",
        dbr.sd=1)
      
      scobj.h.scdbl <- scobj.ds.harmony.sce %>% 
        as.Seurat %>% 
        subset(scDblFinder.class == "singlet")
      
      dim(scobj.ds.harmony)
      # [1] 21885 23113
      dim(scobj.h.scdbl)
      # [1] 21885 23002
    }
    # DoubletFinder---------------------------------------------------------------
    if(FALSE){
      # Identify doublet one sample by one sample
      # https://github.com/chris-mcginnis-ucsf/DoubletFinder
      # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/189
      # https://zhuanlan.zhihu.com/p/469335625
      # https://zhuanlan.zhihu.com/p/668581461
      # https://cloud.tencent.com/developer/article/1825672
      scobj.ds.harmony.split <- SplitObject(scobj.ds.harmony, split.by = "orig.ident") # into list
      
      for (i in 1:length(scobj.ds.harmony.split)) {
        # pK Identification (ground-truth) -------------------------------------------
        sweep.list <- paramSweep(scobj.ds.harmony.split[[i]], PCs = 1:30, sct = TRUE)
        sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
        ## Homotypic Doublet Proportion Estimate -------------------------------------
        homotypic.prop <- modelHomotypic(scobj.ds.harmony.split[[i]]@meta.data$SCT_snn_res.0.4)  
        ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi <- round(0.035 * nrow(scobj.ds.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        ## Run DoubletFinder with varying classification stringencies ----------------
        # scobj.ds.harmony.split[[i]] <- doubletFinder(scobj.ds.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
        scobj.ds.harmony.split[[i]] <- doubletFinder(scobj.ds.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
      }
      
      Singlet <- c(
        rownames(scobj.ds.harmony.split[[1]]@meta.data)[scobj.ds.harmony.split[[1]]@meta.data$DF.classifications_0.25_0.02_112 == "Singlet"],
        rownames(scobj.ds.harmony.split[[2]]@meta.data)[scobj.ds.harmony.split[[2]]@meta.data$DF.classifications_0.25_0.2_108 == "Singlet"],
        rownames(scobj.ds.harmony.split[[3]]@meta.data)[scobj.ds.harmony.split[[3]]@meta.data$DF.classifications_0.25_0.02_112 == "Singlet"],
        rownames(scobj.ds.harmony.split[[4]]@meta.data)[scobj.ds.harmony.split[[4]]@meta.data$DF.classifications_0.25_0.17_106 == "Singlet"],
        rownames(scobj.ds.harmony.split[[5]]@meta.data)[scobj.ds.harmony.split[[5]]@meta.data$DF.classifications_0.25_0.01_110 == "Singlet"])
      
      scobj.ds.harmony@meta.data$id <- rownames(scobj.ds.harmony@meta.data)
      
      scobj.h.dblfinder <- subset(scobj.ds.harmony, subset = id %in% Singlet)
      
      dim(scobj.ds.harmony)
      # [1] 20683 17500
      dim(scobj.h.dblfinder)
      # [1] 20683 16952
    }
  }
  
  scobj@meta.data$barcode = rownames(scobj@meta.data)
  # scobj.ds.sc <- subset(scobj.ds, subset = barcode %in% rownames(scobj.h.scdbl@meta.data))
  scobj.ds.sc <- subset(scobj, subset = barcode %in% rownames(scobj.h.scdbl@meta.data))
  
  # Apply sctransform normalization & Integration & Cluster---------------------
  scobj.ds.harmony <- scobj.ds.sc %>% 
    SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  clustree(scobj.ds.harmony, prefix = "SCT_snn_res.")
  # Visualization-----------------------------------------------------------------
  #UMAP
  scobj.ds.harmony <- RunUMAP(scobj.ds.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
  #T-SNE
  scobj.ds.harmony <- RunTSNE(scobj.ds.harmony, reduction = "harmony", dims = 1:30)
  
  p1 <- DimPlot(scobj.ds.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend()
  p2 <- DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.4", reduction = "umap", label = T)
  p1 + p2
  # Annotation--------------------------------------------------------------------
  scobj.ds.harmony <- RegroupIdents(scobj.ds.harmony, "SCT_snn_res.0.4")

  # curated on the basis of 35657798
  if(TRUE){
    # proximal tubular cells
    PT <- c("MIOX", "ALDOB", "FABP1", "PCK1", "ANPEP") # 0 3
    # tubular progenitor cells (PG)
    PG <- c("CD24", "PROM1", "MMP7") 
    SMC <- c("ACTA2", "RGS5", "MYH11", "TAGLN") # 15 16
    Macrophage <- c("CD68", "CD163", "LYZ") # 18
    Monocyte <- c("CD14", "LYZ", "S100A12", "S100A9", "S100A8")
    DC <- c("FCER1A", "CD1E", "CD1C", "HLA-DMA", "HLA-DMB") 
    Neutrophil <- c('CSF3R', 'S100A8', 'RETNLG') 
    NK <- c("KLRD1", "KLRC1", "GZMB", "PRF1") 
    Fibroblast <- c("SFRP2", "SPARC", "MMP2", "COL3A1", "COL1A1", "COL1A2", 
                    "EMILIN1", "PDGFRB") # 
    EC <- c("PECAM1", "PLVAP", "CDH5", "KDR") # 9 17
    CD8T <- c("CD3D", "CD3E", "CD8A") 
    CD4T <- c("CD3E", "CD3D", "IL7R")  
    T_cell <- c("CD3E", "CD3D", "IL7R", "GZMK")
    B <- c("CD79A", "CD79B", "MS4A1", "CD19") 
    Plasma <- c("IGKC") 
    Mast <- c("TPSAB1", "TPSB2", "KIT")
    CellCycle	<- c('MKI67', 'CCNA2', 'CCNB2', 'PCNA', 'STMN1')
    # RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
  }
  # 37468583
  if(TRUE){
    ####epithelial cells----
    renal_corpuscle <- c('PTPRQ', 'WT1', 'NTNG1', 'NPHS1', 'NPHS2', 'CLIC5', 'PODXL') # 12
    
    glomerulus <- c('CLDN1', 'VCAM1', 'CFH', 'RBFOX1', 'ALDH1A2') # 10
    
    proximal_tubules <- c(
      'LRP2', 'CUBN', 'SLC13A1', 'SLC5A12', 'SLC13A3', 'SLC22A6', 'PRODH2', 'SLC5A2', 
      'SLC22A8', 'SLC34A1', 'SLC22A7', 'MOGAT1', 'SLC5A11', 'SLC22A24', 'SLC7A13', 
      'SLC5A8', 'ABCC3', 'SATB2'
    ) # 0 3
    
    intermediate_tubules <- c(
      'CRYAB', 'TACSTD2', 'SLC44A5', 'KLRG2', 'COL26A1', 'BOC', 'VCAM1', 'SLC39A8', 
      'AQP1', 'LRRC4C', 'LRP2', 'UNC5D', 'SATB2', 'JAG1', 'ADGRL3', 'ID1', 'CLDN1', 
      'AKR1B1', 'CLDN4', 'BCL6', 'SH3GL3', 'SLC14A2', 'SMOC2', 'BCAS1', 'CLCNKA', 
      'CLDN10', 'PROX1'
    ) # 10 ?11
    
    Distal_tubules <- c(
      'CASR', 'SLC12A1', 'UMOD',
      'NELL1', 'ESRRB', 'EGF', 'CLDN14', 'PROX1', 'MFSD4A', 'KCTD16', 'RAP1GAP', 'ANK2', 
      'CYFIP2', 'PPM1E', 'GP2', 'ENOX1', 'TMEM207', 'TMEM52B', 'CLDN16', 'WNK1M',
      'NOS1', 'ROBO2', 'CALCR', 'PPFIA2', 'PAPPA2', 'SLC12A3', 'CNNM2', 'FGF13', 'KLHL3', 
      'LHX1', 'TRPM6', 'TRPM7', 'ADAMTS17', 'ITPKB', 'ZNF385D', 'HS6ST2',
      'TRPV5', 'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1'
    ) %>% unique # 
    # detail--------------------------------------------------------------------
    Thick_Ascending_Limb_Cell <- c(
      'CASR', 'SLC12A1', 'UMOD', 'NELL1', 'ESRRB', 'EGF', 'CLDN14', 'PROX1', 'MFSD4A', 
      'KCTD16', 'RAP1GAP', 'ANK2', 'CYFIP2', 'NELL1', 'ESRRB', 'EGF', 'PPM1E', 'GP2', 
      'ENOX1', 'TMEM207', 'TMEM52B', 'CLDN16', 'WNK1', 'NOS1', 'ROBO2', 'CALCR', 
      'PPFIA2', 'PAPPA2'
    ) %>% unique # ?1 2 ?7 8 ?13 
    Distal_Convoluted_Tubule_Cell <- c(
      'SLC12A3', 'CNNM2', 'FGF13', 'KLHL3', 'LHX1', 'TRPM6', 'TRPM7', 'ADAMTS17', 
      'ITPKB', 'ZNF385D', 'HS6ST2', 'TRPV5', 'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1'
    ) %>% unique # 1 4 13
    
    Collecting_tubules <- c(
      'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1', 'KITLG', 'PCDH7', 'RALYL', 'TOX', 'SGPP1', 
      'SCNN1G', 'SCNN1B', 'KCNIP1', 'GATA3', 'AQP2', 'AQP3', 'FXYD4', 'SOX5', 'PDE10A', 
      'SLC25A29', 'ST6GAL1', 'PAPPA', 'SYK', 'FAM81A', 'PROM1', 'KCNK13', 'FXYD4', 'SOX5',
      'PHACTR1', 'PCDH7', 'SLC14A2', 'HS3ST5', 'TACSTD2', 'TP63', 'GPX2', 'FXYD3', 'KRT5',
      'ATP6V0D2', 'ATP6V1C2', 'TMEM213', 'CLNK', 'SLC4A1', 'SLC26A7', 'HS6ST3', 'NXPH2', 
      'LEF1', 'ADGRF5', 'SLC8A1', 'SCN2A', 'CALB1', 'KIT', 'AQP6', 'STAP1', 'FAM184B', 
      'CALCA', 'SLC4A9', 'SLC35F3', 'SLC26A4', 'INSRR', 'TLDC2'
    ) %>% unique #
    # detail---------------------------------
    Connecting_Tubule_sum <- c('SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1') # 4 5 13
    Connecting_Tubule <- c(
      'KITLG', 'PCDH7', 'RALYL', 'TOX', 
      'SGPP1', 'SCNN1G', 'SCNN1B', 'KCNIP1') %>% unique # 4 5
    Principal_Cell <- c(
      'GATA3', 'AQP2', 'AQP3', 'SCNN1G', 'SCNN1B', 'FXYD4', 'SOX5', 'PDE10A', 
      'SLC25A29', 'ST6GAL1', 'PAPPA', 'SCNN1G', 'SCNN1B', 'FXYD4', 'SOX5', 'SYK', 
      'FAM81A', 'PROM1', 'KCNK13', 'FXYD4', 'SOX5', 'PHACTR1', 'PCDH7', 'SLC14A2', 
      'HS3ST5', 'TACSTD2', 'TP63', 'GPX2', 'FXYD3', 'KRT5') %>% unique # 5
    Intercalated_Cell <- c(
      'ATP6V0D2', 'ATP6V1C2', 'TMEM213', 'CLNK', 'SLC4A1', 'SLC26A7', 'HS6ST3', 
      'NXPH2', 'LEF1', 'ADGRF5', 'SLC4A1', 'SLC26A7', 'SLC8A1', 'SCN2A', 'CALB1', 
      'SLC4A1', 'SLC26A7', 'KIT', 'AQP6', 'STAP1', 'FAM184B', 'CALCA', 'SLC4A9', 
      'SLC35F3', 'SLC26A4', 'INSRR', 'TLDC2') %>% unique # 6 14
    
    #### Endothelial Cell----
    EC <- c('CD34', 'PECAM1', 'PTPRB', 'MEIS2', 'FLT1', 'EMC',"VWF") # 8 16
    #### Vascular Smooth Muscle Cell / Pericyte(stroma cells)----
    VSMC <- c('NOTCH3', 'PDGFRB', 'ITGA8')
    #### Fibroblast(stroma cells)----
    FB <- c('COL1A1', 'COL1A2', 'C7', 'NEGR1', 'FBLN5', 'DCN', 'CDH11') # 13
    #### Immune Cells----
    `B Cell` <- c('BANK1', 'BLK', 'MS4A1', 'BACH2')
    `Plasma Cell` <- c('IGKC', 'TENT5C', 'MZB1', 'FCRL5', 'CD38', 'JCHAIN')
    `T Cell` <- c('CD96', 'CD247', 'THEMIS', 'BCL11B', 'CAMK4', 'IL7R')
    `Natural Killer T Cell` <- c(
      'CD96', 'CD247', 'RUNX3', 'GNLY', 'NKG7', 'CCL5', 'KLRF1', 'CCL4', 'GZMA')
    `Mast Cell` <- c('MS4A2', 'CPA3', 'KIT')
    `M2 Macrophage` <- c('F13A1', 'MRC1', 'CD163', 'STAB1', 'SLC1A3', 'CD14', 'FOLR2')
    `Monocyte-derived Cell` <- c(
      'MSR1', 'ITGAX', 'HLA-DQA1', 'HLA_DRB1', 'CSF2RA', 'CD14', 'TRPM2') # ?17
    `Classical Dendritic Cell` <- c(
      'ITGAX', 'HLA-DQA1', 'HLA-DRA', 'CSF2RA','CIITA', 'WDFY4', 'FLT3', 
      'ZNF366', 'CADM1', 'ZBTB46', 'CLEC9A') # 17
    `Plasmacytoid Dendritic Cell` <- c('IRF8', 'CUX2', 'P2RY14', 'IL3RA', 'CLEC4C')
    `Non-Classical Monocyte` <- c(
      'CTSS', 'IRAK3', 'TCF7L2', 'TNFRSF1B', 'FCN1', 'HLA-DRA', 'FCGR3A') # 17
    Neutrophil <- c('S100A9', 'S100A8', 'IFITM2', 'FCGR3B', 'CD1C')
  }
  
  VlnPlot(scobj.ds.harmony, features = Thick_Ascending_Limb_Cell, group.by = "SCT_snn_res.0.4", assay = "RNA", layer = "counts", log = TRUE)
  
  scobj.ds.harmony@meta.data$cell_type <- "Unknown"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(0,3)] <- "proximal_tubules"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(1,4,13)] <- "Distal_Convoluted_Tubule_Cell"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(2,7,8)] <- "Thick_Ascending_Limb_Cell"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(4,5)] <- "Connecting_Tubule & Principal_Cell"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(6,14)] <- "Intercalated_Cell"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(9,17)] <- "EC"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(10,11)] <- "glomerulus & intermediate_tubules"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(12)] <- "renal_corpuscle"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(15,16)] <- "SMC"
  scobj.ds.harmony@meta.data$cell_type[scobj.ds.harmony@meta.data$SCT_snn_res.0.4 %in% c(18)] <- "Macrophage"
  
  p1 <- DimPlot(scobj.ds.harmony, reduction = "umap",
                group.by = "SCT_snn_res.0.4",
                label = TRUE, pt.size = 0.5) 
  p2 <- DimPlot(scobj.ds.harmony, reduction = "umap",
                group.by = "cell_type",
                label = TRUE, pt.size = 0.5)
  p3 <- DimPlot(scobj.ds.harmony, reduction = "umap",
                group.by = "orig.ident",
                label = TRUE, pt.size = 0.5)+ NoLegend()
  # p4 <- FeaturePlot(scobj.ds.harmony, reduction = "umap", features = "FGL1")
  p1 + p2 + p3 + plot_layout(nrow = 2, ncol = 2)
  }
## snATACseq--------------------------------------------------------------------
if(TRUE){
  # Subsequently datasets were processed with Seurat v3.0.2 and its companion package 
  # Signac v0.2.1 (https://github.com/timoast/signac)21. 
  
  # Low-quality cells were removed from the aggregated snATAC-seq library
  # (peak region fragments > 2500, peak region fragments < 25000, %reads in peaks > 15, 
  # blacklist ratio < 0.001, nucleosome signal < 4 & mitochondrial gene ratio < 0.25) 
  
  # before normalization with term-frequency inverse-document-frequency (TFIDF).
  
  # A fraction of reads in peaks, number of reads in peaks per cell and ratio reads 
  # in genomic blacklist region per cell for each patient were shown in Supplementary
  # Fig. 19. 
  
  # Dimensional reduction was performed via singular value decomposition (SVD) of 
  # the TFIDF matrix and UMAP. 
  
  # A KNN graph was constructed to cluster cells with the Louvain algorithm. 
  
  # Batch effect was corrected with Harmony67 using the “RunHarmony” function in Seurat. 
  
  # A gene activity matrix was constructed by counting ATAC peaks within the gene 
  # body and 2 kb upstream of the transcriptional start site using protein-coding 
  # genes annotated in the Ensembl database. 
  
  # The gene activity matrix was log-normalized prior to label
  # transfer with the aggregated snRNA-seq Seurat object using canonical correlation
  # analysis. 
  
  # The aggregated snATAC-seq object was filtered using a 97% confidence
  # threshold for cell-type assignment following label transfer to remove heterotypic
  # doublets. 
  
  # The filtered snATAC-seq object was reprocessed with TFIDF, SVD, and
  # batch effect correction followed by clustering and annotation based on lineage 
  # specific gene activity. After filtering, there was a mean of 5408 ± 1393 nuclei per
  # snATAC-seq library with a mean of 7538 ± 2938 peaks detected per nucleus. The
  # final snATAC-seq library contained a total of 214,890 unique peak regions among
  # 27,034 nuclei and represented all major cell types within the kidney cortex 
  # (Supplementary Table 2, Supplementary Fig. 20). 
  
  # Differential chromatin accessibility between cell types was assessed with the 
  # Signac FindMarkers function for peaks detected in at least 20% of cells using a 
  # likelihood ratio test and a log-fold-change threshold of 0.25. Bonferroni-adjusted 
  # p-values were used to determine significance at an FDR < 0.05. 
  
  # Genomic regions containing snATAC-seq peaks were annotated
  # with ChIPSeeker68 (v1.24.0) and clusterProfiler69 (v3.16.1) using the UCSC
  # database70 on hg38.
}














