# Infercnv======================================================================
# https://github.com/broadinstitute/inferCNV/wiki/File-Definitions
# https://github.com/broadinstitute/inferCNV/wiki/instructions-create-genome-position-file
# https://github.com/broadinstitute/infercnv/wiki/Output-Files
# https://github.com/broadinstitute/infercnv/issues/231
# https://github.com/broadinstitute/infercnv/issues/520
# https://bioc.r-universe.dev/infercnv/doc/manual.html

# https://blog.csdn.net/zfyyzhys/article/details/139630834
# https://statbiomed.github.io/SingleCell-Workshop-2021/CNV-analysis.html
# There are four major methods for CNV analyis in scRNA-seq:
#   
# inferCNV
# CopyKat by Gao et al 2020 #https://github.com/navinlabcode/copykat
# CaSpER by Harmanci et al, 2020 
# https://www.bioconductor.org/packages/release/bioc/html/casper.html
# https://github.com/akdess/CaSpER
# HoneyBadger by Fan et al

if(infercnv){
  # inferCNV object files
  # During the inferCNV data processing, it will write persistent object files representing the inferCNV object state at the corresponding stage of the processing. Each processing stage is numbered, so it's straightforward to follow the operations performed.
  # These objects can be resurrected in an R session like so:
  # library(infercnv)
  # infercnv_obj = readRDS('preliminary.infercnv_obj')
  # The above 'preliminary' object is used to generate the 'infercnv.preliminary.png' file, and the data included are the target for either denoising or performing CNV prediction using HMMs.
  # If you want to interact with any of the object files, the essential features to be aware of include the following slots of the inferCNV S4 object:
  
  # 'infercnv_obj@ expr.data' : contains the processed expression matrix as it exists at the end of that stage for which that inferCNV object represents.
  # 
  # 'infercnv_obj@reference_grouped_cell_indices' : list containing the expression matrix column indices that correspond to each of the normal (reference) cell types.
  # 
  # 'infercnv_obj@observation_grouped_cell_indices' : similar list as above, but corresponds to the tumor cell types.
  # 
  # Based on the above slots, it would be straightforward to extract info of interest and/or move data into other analysis frameworks.
  
  # https://mp.weixin.qq.com/s/RFHyADX1mtpDGKNE7hDEbg
  # https://www.jianshu.com/p/815d8acb2ada
  # https://zhuanlan.zhihu.com/p/392282792
  
  # Step 1----------------------------------------------------------------------
  # run on Linux
  # refdata-gex-GRCh38-2024-A
  python gtf_to_position_file.py --attribute_name gene_name genes.gtf your_gen_pos.txt
  
  library(dplyr)
  library(Seurat)
  library(infercnv)
  
  setwd("/home/User/")
  set.seed(1011)
  options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB
  
  load("scobj.harmony.Rdata")
  
  scobj.harmony <- subset(scobj.harmony, 
                          subset = cell_type %in% c("B cell","T cell", "Epithelial"))
  
  counts <- scobj.harmony@assays$RNA$counts
  
  # a description of the cells, indicating the cell type classifications
  # Sample annotation file
  # The sample annotation file is used to define the different cell types, and optionally, indicating how the cells should be grouped according to sample (ie. patient). The format is simply two columns, tab-delimited, and there is no column header.
  # 
  # MGH54_P2_C12    Microglia/Macrophage
  # MGH36_P6_F03    Microglia/Macrophage
  # MGH54_P16_F12   Oligodendrocytes (non-malignant)
  # MGH54_P12_C10   Oligodendrocytes (non-malignant)
  # MGH36_P1_B02    malignant_MGH36
  # MGH36_P1_H10    malignant_MGH36
  # The first column is the cell name, and the 2nd column indicates the known cell type. For the normal cells, if you have different types of known normal cells (ie. immune cells, normal fibroblasts, etc.), you can give an indication as to what the cell type is. Otherwise, you can group them all as 'normal'. If multiple 'normal' types are defined separately, the the expression distribution for normal cells will be explored according to each normal cell grouping, as opposed treating them all as a single normal group. They'll also be clustered and plotted in the heatmap according to normal cell grouping.
  # The sample (ie. patient) information is encoded in the attribute name as "malignant_{patient}", which allows the tumor cells to be clustered and plotted according to sample (patient) in the heatmap.
  annot  <- scobj.harmony@meta.data %>%
    mutate(cell_type = case_when(
      cell_type == "Epithelial" ~ paste0(sample,"_",cell_type),
      TRUE ~ cell_type)) %>%
    dplyr::select(cell_type)
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix= counts,
    gene_order_file="your_gen_pos.txt",
    annotations_file = annot,
    delim="\t",
    ref_group_names=c("B cell","T cell")) 
  
  # Please use "options(scipen = 100)" before running infercnv if you are using the analysis_mode="subclusters" option or you may encounter an error while the hclust is being generated.
  options(scipen = 100)
  
  infercnv_obj_run = infercnv::run(
    infercnv_obj, # An infercnv object populated with raw count data
    cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir = "/home/User/infercnv/",
    analysis_mode='subclusters',
    # tumor_subcluster_partition_method = "random_trees", # refer to 37832554 # very very slow
    num_threads = 12,
    # leiden_resolution = 0.1,
    cluster_by_groups = TRUE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
    output_format="pdf",
    write_expr_matrix = TRUE, # important
    num_threads = 8,
    denoise = TRUE, # denoise，sd_amplifier，noise_logistic
    HMM = FALSE # very very slow
    )
  # save(infercnv_obj_run, file = "infercnv.Rdata") 
  
  # nohup Rscript infercnv.R &
  # Step 2----------------------------------------------------------------------
  # Ref
  if(FALSE){
    # 2025, NG, pmid: 40169792
    # Spatial CNV analysis
    # The Python-based tool infercnvpy (https://github.com/icbi-lab/infercnvpy)
    # was used to deduce the copy number variation (CNV) status in
    # each spot. In brief, 
    # spots that were histologically normal and identified as type 1 luminal cells 
    # through spatial transcriptomics expression analysis were used as references 
    # for CNV computations. 
    # To avoid batch discrepancies, the references used for each patient-derived 
    # sample  were sourced from the corresponding individual. 
    # Spots that were confirmed as stromal in their expression profile and 
    # morphologically validated were excluded from the subclonal analysis. 
    # Upon obtaining the CNV score matrix, we sequentially applied cnv.tl.pca, cnv.
    # pp.neighbors and cnv.tl.leiden to discern different tumor subclones.
    # This approach is aligned with methods that have been verified as effective
    # in previous studies31,65.
    
    # Note, the current inferCNV software generates a lot of output files, and only a small subset are of greatest interest. In the future, outputs will be better organized to facilitate identification of the most useful outputs.
    # infercnv.preliminary.png : the preliminary inferCNV view (prior to denoising or HMM prediction)
    # infercnv.png : the final heatmap generated by inferCNV with denoising methods applied.
    # infercnv.references.txt : the 'normal' cell matrix data values.
    # infercnv.observations.txt : the tumor cell matrix data values
    # infercnv.observation_groupings.txt : group memberships for the tumor cells as clustered.
    # infercnv.observations_dendrogram.txt : the newick formatted dendrogram for the tumor cells that matches the heatmap.
  
    # If you want to interact with any of the object files, the essential features to be aware of include the following slots of the inferCNV S4 object:
    #   
    # 'infercnv_obj@ expr.data' : contains the processed expression matrix as it exists at the end of that stage for which that inferCNV object represents.
    # 'infercnv_obj@reference_grouped_cell_indices' : list containing the expression matrix column indices that correspond to each of the normal (reference) cell types.
    # 'infercnv_obj@observation_grouped_cell_indices' : similar list as above, but corresponds to the tumor cell types.
    # 
    # Based on the above slots, it would be straightforward to extract info of interest and/or move data into other analysis frameworks.
    }
  # infercnv.references.txt : the 'normal' cell matrix data values.
  infercnv.references <- data.table::fread("G:/infercnv.references.txt")
  # infercnv.observations.txt : the tumor cell matrix data values
  infercnv.observation <- data.table::fread("G:/infercnv.observations.txt")
  # infercnv.observation_groupings.txt : group memberships for the tumor cells as clustered.
  infercnv.observation_groupings <- data.table::fread("G:/infercnv.observation_groupings.txt")
  
  preliminary = readRDS('G:/preliminary.infercnv_obj')
  run.final = readRDS('G:/run.final.infercnv_obj')
  
  dim(infercnv.references)
  # [1]  5137 42613
  head(infercnv.references[,1:6])
  
  dim(infercnv.observation)
  # [1] 5137 8392
  head(infercnv.observation[,1:6])
  
  # √ Method 1 -------------------------------------------------------------------
  # 2025, Nat Aging, pmid: 40211000
  # ..., we integrated the results of spike-in (an example with EOPC1 shown in Fig. 1f) 
  # and copy number variation (CNV) scores based on inferCNV, which is an algorithm 
  # inferring CNV levels with transcriptomic matrices, and classified epithelial 
  # cells into three subtypes: malignant, normal, and other epithelia.
  # 
  # Malignant epithelial cells identification. 
  # Based on the CNV of epithelial
  # cells inferred by inferCNV (v.1.10.1) package with immune cells
  # as reference, we identified malignant epithelia through two methods:
  #   spike-in and CNV scores. 
  # For spike-in method, we mixed test cells with reference cells and 
  # conducted k-means clustering based on inferred CNV levels; 
  # test cells that cluster together with reference cells represent
  # normal epithelial cells. 
  # CNV score of each cell cluster was calculated with 
  # Σ (CNVi − 1)2 , 
  # where CNVi represents the CNV level inferred for region i.
  CNV_score <- function(obj = infercnv.observation){
    obj <- obj[,-1]
    res <- colSums((obj - 1)^2)
    return(res)
  }
  # √ Method 2 -------------------------------------------------------------------
  # 2025, Cancer Discov, pmid: 39774838
  # we inferred sCNAs from scRNA-seq tumor epithelial cells from 29 tumor samples 
  # using CopyKAT and validated the sCNA calls using inferCNV (median rho = 0.55 
  # using Spearman’s rank correlation, P < 2.20 × 10−16 in 29 samples).
  
  # Method 3 (?)----------------------------------------------------------------
  # 2025, Exp Mol Med, pmid: 39741182
  # CNV estimation
  # We used R package infercnv (v0.8.2) to identify somatic copy number
  # variations and we used epithelial cells of adjacent tissues as the reference.
  # Each cell for the extent of CNV signal was scored and the mean of squares
  # of CNV values was defined across the genome. We then defined the
  # putative malignant cells as those with CNV signal more than 0.05 and CNV
  # correlation more than 0.5.
  # CNV_score3 <- function(obj = infercnv.observation){
  #   obj <- obj[,-1]
  #   res <- colMeans(obj^2)
  #   return(res)
  # }
  # Method 4 (?)----------------------------------------------------------------
  # 2023, Cancer Cell, pmid: 37832554
  # Single-cell copy number variation analysis
  # Copy number instability was assessed with the R package inferCNV (version 1.9.1), 
  # which is designed to infer copy number variation (CNV) from single-cell RNA-seq 
  # data. This package compares the expression intensities of genes across epithelial
  # cells and relates this to expression in normal cells. A random subset of T cells 
  # and B cells (5000 of each) was set as the reference normal cells, while
  # sex chromosomes were excluded as described in the references to this manuscript. 
  # All epithelial cells were clustered by the tree based algorithm from the 
  # imputed CNV pattern, and the result was visualized on heatmap. Normal epithelial 
  # cells were determined by extremely lower CNV levels and lower expression of 
  # ESCC cancer-related genes (EPCAM, TSTA3, NECTIN4, SFN, KRT5, KRT16, KRT18, 
  # KRT17, KRT19, KRT6A, and KRT6B),67 compared with malignant cells (Figures S1E 
  # and S1F).
  # Method 5 without analysis --------------------------------------------------
  # 2023, Adv Sci, pmid: 36709495
  # 2.5. Copy Number Variants (CNV) Estimation
  # To explore the heterogeneity of copy number alterations existing in epithelial 
  # cells of primary tumor and lymph node metastasis site of ESCC, inferCNV (
  #   https://github.com/broadinstitute/inferCNV) was applied to estimate the changes 
  # of somatic alterations of large-scale chromosomal copy number variants
  # (amplifications or deletions) in a single epithelial cell. 
  # The epithelial cell raw single-cell gene expression data were extracted
  # from the Seurat object followed by the software manual. 
  # Single cell data from epithelial cells in normal tissue were used as
  # a control reference and then the inferCNV analysis with the
  # default parameters was performed.
  
  # Compared to the reference data of normal epithelial cells
  # (ESCC_N: C11, C18), the primary and lymph node metastatic
  # malignant epithelial cells presented a large-scale diverse in chromosomal
  # copy number alteration patterns. For example, amplification
  # of chromosome 1 and deletion of 19 in specific regions
  # could be frequently detected in both ESCC_PT and ESCC_LNM,
  # whiles amplifications of chromosomes 7 and 8 were only observed
  # in ESCC_PT and ESCC_LNM, respectively (Figure 5b).
  # √ Method 6 -------------------------------------------------------------------
  # 2022, Cancer Lett, pmid: 35917973
  # 2.7. InferCNV analysis
  # We isolated ductal cells and construct a new gene-cell matrix. Somatic
  # large-scale chromosomal CNV score of each ductal cell was
  # calculated using the R package inferCNV (v1.6.0). A raw counts matrix,
  # annotation file, and gene/chromosome position file were prepared according
  # to data requirements (https://github.com/broadinstitute/inferCNV). 
  # T cells, macrophages, endothelial cells, and stellate cells were
  # selected as reference normal cells. 
  # The default parameters were applied (cutoff = 0; denoise = 0.2). 
  # The CNV score was calculated as quadratic sum of CNV region.
  CNV_score6 <- function(obj = infercnv.observation){
    obj <- obj[,-1]
    res <- colSums(obj^2)
    return(res)
  }
  # Method 7 without analysis --------------------------------------------------
  # inferCNV analysis. 
  # To identify malignant cells, we identified evidence for somatic
  # alterations of large-scale chromosomal copy number variants, either gains or losses,
  # in a single cell using inferCNV (https://github.com/broadinstitute/inferCNV), in
  # addition to the expression of EPCAM. The raw single-cell gene expression data was
  # extracted from the Seurat object according to the software recommendation. 
  # A public single-cell data derived from normal epithelium cells was included as a
  # control reference32 (GEO accession number: GSE121600). 
  # We preformed inferCNV analysis with the default parameters.
  
  
}
if(CopyKAT){
  # https://www.jianshu.com/p/78089451cc0e
  # https://www.jianshu.com/p/4274ea0a9eec
  # 最后作者提到一个需要注意的点，不是所有的肿瘤都存在CNV。儿童肿瘤和血液肿瘤中基本没有copy number event，所以是不适合用这些方法（copyKAT或inferCNV）来寻找肿瘤细胞的。
  # https://cloud.tencent.com/developer/article/1910795
  # https://blog.csdn.net/weixin_53637133/article/details/138145982
  library(Seurat)
  library(copykat)
  
  tmp <- subset(scobj.harmony, subset = cell_type == "Epithelial")
  exp.rawdata <- as.matrix(tmp@assays$RNA$counts)
  # write.table(exp.rawdata, file="exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)
  
  copykat.test <- copykat(
    rawmat=exp.rawdata, 
    id.type="S", 
    ngene.chr=5, 
    win.size=25, 
    KS.cut=0.1, # 增加KS.cut会降低敏感度，通常范围在0.05-0.15
    sam.name="test", 
    distance="euclidean", 
    norm.cell.names="",
    output.seg="FLASE", 
    plot.genes="TRUE", 
    genome="hg20",
    n.cores=1)
  
  table(copy$prediction_0.1, copy$prediction_0.05)
  # aneuploid diploid not.defined
  # aneuploid         968       2           0
  # diploid             1    6541           0
  # not.defined         0       0        1528
  
  prediction <- data.table::fread("test_copykat_prediction.txt")
  
  scobj.harmony@meta.data$cell.names <- rownames(scobj.harmony@meta.data)
  meta.data <- scobj.harmony@meta.data %>% 
    left_join(prediction, by = "cell.names")
  rownames(meta.data) <- meta.data$cell.names
  scobj.harmony@meta.data <- meta.data
  
  
  DotPlot(subset(scobj.harmony, subset = cell.names %in% prediction$cell.names),
          features = c("EPCAM","TSTA3","NECTIN4","SFN","KRT5","KRT16","KRT18","KRT17","KRT19","KRT6A","KRT6B"), 
          group.by = "copykat.pred") + 
    scale_color_viridis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
