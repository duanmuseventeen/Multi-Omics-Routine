rm(list = ls());gc()

# Load pkgs---------------------------------------------------------------------
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)
# library(decontX)
library(scDblFinder)
library(DoubletFinder)
library(Augur)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(patchwork)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(GSVA)
# library(slingshot) # Trajectory
# library(tradeSeq) # Trajectory
# library(monocle) # Trajectory
library(monocle3) # Trajectory
library("CellChat") # Communication
# library("iTALK") # Communication
# library(infercnv) # CNV
# library(copykat) # CNV
# Load data---------------------------------------------------------------------
set.seed(1011)
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')


pdac_a <- data.table::fread("GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz")
pdac_b <- data.table::fread("GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz")
# barcodes <- data.table::fread("GSE111672_1000L3_barcodes.txt.gz")
# barcode1 <- data.table::fread("GSE111672_inDrop_gel_barcode1_list.txt.gz")    
# barcode2 <- data.table::fread("GSE111672_inDrop_gel_barcode2_list.txt.gz")     
GSM3036909 <- data.table::fread("GSM3036909_PDAC-A-indrop1.tsv.gz")              
GSM3036910 <- data.table::fread("GSM3036910_PDAC-A-indrop2.tsv.gz")              
GSM3405527 <- data.table::fread("GSM3405527_PDAC-A-indrop3.tsv.gz")              
GSM3405528 <- data.table::fread("GSM3405528_PDAC-A-indrop4.tsv.gz")              
GSM3405529 <- data.table::fread("GSM3405529_PDAC-A-indrop5.tsv.gz")              
GSM3405530 <- data.table::fread("GSM3405530_PDAC-A-indrop6.tsv.gz")              
GSM3405531 <- data.table::fread("GSM3405531_PDAC-B-indrop1.tsv.gz")              
GSM3405532 <- data.table::fread("GSM3405532_PDAC-B-indrop2.tsv.gz")              
GSM3405533 <- data.table::fread("GSM3405533_PDAC-B-indrop3.tsv.gz")              
GSM4100717 <- data.table::fread("GSM4100717_PDAC-C-indrop1.tsv.gz")              
GSM4100718 <- data.table::fread("GSM4100718_PDAC-C-indrop2.tsv.gz")              
GSM4100719 <- data.table::fread("GSM4100719_PDAC-C-indrop3.tsv.gz")              
GSM4100720 <- data.table::fread("GSM4100720_PDAC-C-indrop4.tsv.gz") 


table(colnames(pdac_a))
# Acinar cells                         Cancer clone A 
# 13                                    126 
# Cancer clone B            Ductal - APOL1 high/hypoxic 
# 170                                    215 
# Ductal - CRISP3 high/centroacinar like                  Ductal - MHC Class II 
# 529                                    287 
# Ductal - terminal ductal like                        Endocrine cells 
# 350                                      3 
# Endothelial cells                            Fibroblasts 
# 11                                      5 
# Genes                          Macrophages A 
# 1                                     21 
# Macrophages B                             Mast cells 
# 19                                     14 
# mDCs A                                 mDCs B 
# 12                                     33 
# Monocytes                                   pDCs 
# 18                                     13 
# RBCs                     T cells & NK cells 
# 15                                     40 
# Tuft cells 
# 32 

table(colnames(pdac_b))
# Acinar cells                         Cancer clone A 
# 6                                    339 
# Ductal - CRISP3 high/centroacinar like                  Ductal - MHC Class II 
# 152                                    211 
# Ductal - terminal ductal like                        Endocrine cells 
# 736                                     13 
# Endothelial cells                                  Genes 
# 159                                      1 
# Macrophages                             Mast cells 
# 9                                     13 
# mDCs                              Monocytes 
# 35                                     20 
# RBCs                             Tuft cells 
# 3                                     37 






