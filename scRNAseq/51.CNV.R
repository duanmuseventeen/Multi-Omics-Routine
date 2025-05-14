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