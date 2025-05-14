# Ref:
# https://github.com/aertslab/pySCENIC
# https://pyscenic.readthedocs.io/en/latest/
# https://scenic.aertslab.org/
# http://scope.aertslab.org

# Input: expression matrix
# The input to SCENIC is the single-cell RNA-seq expression matrix:
# Each column corresponds to a sample (cell) and each row corresponds to a gene.
# The gene ID should be the gene-symbol and stored as rownames (for compatibility with RcisTarget annotation databases).
# Expression units: The preferred expression values are gene-summarized counts. There is currently not a strong recommendation towards using the raw counts, or counts normalized through single-cell specific methods (e.g. Seurat). Other measurements, such as transcripts/counts per million (TPM) and FPKM/RPKM, are also accepted as input. However, note that some authors recommend avoiding within sample normalization (i.e. TPM) for co-expression analysis (first step of SCENIC) because they may induce artificial co-variation (Crow et al. (2016)). The choice of input expression matrix might have some effect on the co-expression analysis to create the regulons (step 1). The other steps of the workflow are not directly affected by the input expression values: (2) The expression is not taken into account for the motif analysis, and (3) AUCell, which is used for scoring the regulons on the cells, is cell ranking-based (it works as an implicit normalization). Overall, SCENIC is quite robust to this choice, we have applied SCENIC to datasets using raw (logged) UMI counts, normalized UMI counts, and TPM and they all provided reliable results (see Aibar et al. (2017)).

# Species-specific databases
# In addition to the R-packages, you will also need to download the species-specific databases for RcisTarget (the motif rankings). The links to all the available databases are available in our website. By default, SCENIC uses the databases that score the motifs in the promoter of the genes (up to 500bp upstream the TSS), and in the 20kb around the TSS (+/-10kbp).

# Sample dataset
# Running SCENIC in a real dataset typically takes a few hours. The toy example used in these tutorials (200 cells and <1000 genes) is a subset of a dataset of 3005 cells from the adult mouse brain, including neurons (e.g. pyramidal neurons and interneurons) and glia (oligodendrocytes, astrocytes/ependymal, endothelial/mural and microglia). The expression values are Unique Molecular Identifier counts.
# Zeisel, A., et al. (2015). Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142. doi: 10.1126/science.aaa1934
# The output files from the run on the full dataset are available at http://scenic.aertslab.org/examples/

# RcisTarget数据库文件下载链接：
# 
# TF：https://github.com/aertslab/pySCENIC/tree/master/resources
# Feather：motif-gene rank，https://resources.aertslab.org/cistarget/databases
# Tbl： motif-TF annotation，https://resources.aertslab.org/cistarget/motif2tf
# 
# 作者：生信云笔记
# 链接：https://www.jianshu.com/p/45958b75ec36
# 来源：简书
# 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。

# 2025, NG, pmid: 40169792
# ... The output of the pySCENIC workflow was then transferred into the R environment. Differentially
# activated regulons in each subpopulation were identified by a linear model from the LIMMA package58. 
# The P value was adjusted by the Benjamini–Hochberg procedure. Regulons with an adjusted P
# value < 0.001 were deemed significantly differentially activated.

library(SCENIC)
library(Seurat)
library(ggplot2)
library(dplyr)
# ================================== R =========================================
# ================================== Input =====================================
# Expression matrix
# 1772066100_D04 1772063062_G01 1772060224_F07 1772071035_G09 1772067066_E12
# Arhgap18              0              0              0              0              0
# Apln                  0              0              0              0              0
# Cnn3                  0              0              0              0              0
# Eya1                  0              0              0              4              0
# Mtss1l                0              0              0              0              0
# Boc                   0              0              0              0              0

# Cell info/phenodata
##                        cell_type nGene nUMI
## 1772066100_D04     interneurons   170  509
## 1772063062_G01 oligodendrocytes   152  443
## 1772060224_F07        microglia   218  737
## 1772071035_G09    pyramidal CA1   265 1068
## 1772067066_E12 oligodendrocytes    81  273
## 1772066100_B01    pyramidal CA1   108  191

scenic.obj # a Seurat object
exprMat <- scenic.obj@assays$RNA$counts
cellInfo <- scenic.obj@meta.data
cellInfo$nGene <- colSums(exprMat > 0)
# saveRDS(cellInfo, file="G:/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(cell_type=c("cDC1_CCR7"="forestgreen", 
                           "cDC2_FCER1A"="darkorange", 
                           "pDC_LILRA4"="magenta4"))
colVars$cell_type <- colVars$cell_type[intersect(names(colVars$cell_type), 
                                               cellInfo$cell_type)]
# saveRDS(colVars, file="G:/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$cell_type, legend=names(colVars$cell_type))
# ========================== Initialize SCENIC settings ========================
# In order to keep consistent settings across the multiple steps of SCENIC, most functions in SCENIC package use a common object where the options for the current run are stored. This object replaces the “arguments” for most functions, and should be created at the begining of a SCENIC run with the function initializeScenic().
# https://github.com/aertslab/SCENIC/issues/168
# library(arrow)
# install_arrow(verbose=TRUE)

# try
# trace(checkAnnots, edit = T)
# and edit the fifth line from
# rnktype = "features"
# to
# rnktype = "motifs"

org <- "hgnc" # or mgi, or dmel
# db_mcVersion <- 'v9'
dbDir <- "G:/Database/TF-SCENIC + RcisTarget/human" # RcisTarget databases location
myDatasetTitle <- "SCENIC" # choose a name for your analysis
data(defaultDbNames)
# dbs <- defaultDbNames[[org]]
dbs <- c("hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
         "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")
names(dbs) <- c("500bp","10kb")

# # https://github.com/aertslab/SCENIC/issues/364
# #load in the motifannotation this will load it into your environment but as the name in which is given to the list argument
# data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
# #rename the motif annnotion by attributing it to the variable that is in the error
# motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(
  org=org, 
  dbDir=dbDir, 
  dbs=dbs, 
  datasetTitle=myDatasetTitle, 
  nCores=2) 

motifAnnotations_hgnc <- motifAnnotations

scenicOptions <- initializeScenic(
  org=org, 
  dbDir=dbDir, 
  dbs=dbs, 
  datasetTitle=myDatasetTitle, 
  nCores=2) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "G:/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "G:/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="G:/scenicOptions.Rds") 

...
# ================================== Python ====================================
load("DC.Rdata")
write.csv(DC@assays$RNA$counts, "DC.csv")

# SCENIC steps
# STEP 1: Gene regulatory network inference, and generation of co-expression modules
# Phase Ia: GRN inference using the GRNBoost2 algorithm
# For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system. We use the counts matrix (without log transformation or further processing) from the loom file we wrote earlier. Output: List of adjacencies between a TF and its targets stored in ADJACENCIES_FNAME.

# transcription factors list
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_hg38.txt" # human
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_dmel.txt" # drosophila
# f_tfs = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/resources/allTFs_mm.txt"   # mouse
# tf_names = load_tf_names( f_tfs )

# Step 1
pyscenic grn DC.csv allTFs_hg38.txt -o adj.csv --num_workers 8
# Step 2-3
pyscenic ctx adj.tsv \
hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname DC.csv \
--output reg.csv \
--mask_dropouts \
--num_workers 8
# Step 4
# It is important to check that most cells have a substantial fraction of expressed/detected genes in the calculation of the AUC. The following histogram gives an idea of the distribution and allows selection of an appropriate threshold. In this plot, a few thresholds are highlighted, with the number of genes selected shown in red text and the corresponding percentile in parentheses). See the relevant section in the R tutorial for more information.
# 
# By using the default setting for --auc_threshold of 0.05, we see that 1192 genes are selected for the rankings based on the plot below.
# 
# nGenesDetectedPerCell = np.sum(adata.X>0, axis=1)
# percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
# print(percentiles)
# 0.01     473.58
# 0.05    1192.00
# 0.10    1390.90
# 0.50    1939.00
# 1.00    3998.00
# dtype: float64
# fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
# sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
# for i,x in enumerate(percentiles):
#   fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
# ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small',rotation_mode='anchor' )
# ax.set_xlabel('# of genes')
# ax.set_ylabel('# of cells')
# fig.tight_layout()

pyscenic aucell \
DC.csv \
reg.csv \
--output pyscenic_DC_res.loom pyscenic_DC_res.csv \
--num_workers 8

# Q ----------------------------------------------------------------------------
# 1. pySCENIC运行的原理











