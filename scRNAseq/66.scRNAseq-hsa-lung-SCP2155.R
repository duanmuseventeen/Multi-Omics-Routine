# 2023
# Sci Immunol
# 37027478
# Identification of a broadly fibrogenic macrophage subset induced by type 3 inflammation

# COPD  -   chronic obstructive pulmonary disease
# IPF   -   idiopathic pulmonary fibrosis 
# SSc   -   Systemic Sclerosis

# ATI/II-   type I/II pneumocyte

rm(list = ls());gc()
setwd("~/SCP2155/expression/")
set.seed(1011)
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
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
# Load data---------------------------------------------------------------------
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB
colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

scobj <- "SCP2155" %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "SCP2155",
    assay = "RNA")

meta.data = data.table::fread("~/SCP2155/2022_SI_human_lung_allcells_metadata_file_updated.txt")
meta.data <- as.data.frame(meta.data)
rownames(meta.data) <- meta.data$NAME

all(rownames(scobj@meta.data) == rownames(meta.data))

scobj@meta.data <- meta.data
save(scobj, file = "SCP2155_processed.Rdata")

# DEA --------------------------------------------------------------------------
# cell - cell ==================================================================
myscobj <- subset(scobj, subset = disease__ontology_label %in% c("chronic obstructive pulmonary disease", "normal"))
myscobj <- subset(myscobj, subset = cell_type__ontology_label != "apoptosis fated cell")
myscobj <- subset(myscobj, subset = cell_type__ontology_label != "mitotic cell cycle")

# table(myscobj$health)
# COPD healthy 
# 68792  259912 

myscobj@assays$RNA$data = myscobj@assays$RNA$counts
all.markers <- FindAllMarkers(myscobj, group.by = "cell_type__ontology_label")

save(all.markers, file = "findallmarkers.Rdata")
write.csv(all.markers, file = "findallmarkers.csv")
# health - COPD ================================================================
suppressWarnings(require(crayon, warn.conflicts = F))

Pa    = crayon::cyan
Er    = crayon::red$bold
Sa    = crayon::blue
No    = crayon::magenta$bold
Wa    = crayon::yellow

mymeta <- myscobj@meta.data %>% 
  dplyr::select(biosample_id, donor_id, health) %>% 
  distinct(biosample_id, donor_id, health, .keep_all = TRUE) %>% 
  mutate(biosample_id = str_replace_all(biosample_id, "_", "-"))

for (i in unique(myscobj$cell_type__ontology_label)) {
  tmp <- myscobj %>% subset(subset = cell_type__ontology_label == i)
  
  counts   <- AggregateExpression(tmp, group.by = "biosample_id")
  counts   <- counts$RNA %>% as.data.frame
  colnames(counts) <- str_remove(colnames(counts), "^g")
  mymetacon<- mymeta %>% filter(biosample_id %in% colnames(counts))
  counts   <- counts[,mymetacon$biosample_id]
  condition<- factor(mymetacon$health, levels = c("healthy","COPD"))
  
  message(Sa(paste0("--> Start ", i," DEA <--")))
  if(all(mymetacon$biosample_id == colnames(counts))){
    DEG <- data.frame(
      symbol = rownames(counts),
      health = rowMeans(counts[,condition == "healthy"]),
      copd   = rowMeans(counts[,condition == "COPD"]),
      t      = apply(counts, 1, function(x){t.test(x ~ condition, var.equal = TRUE)$p.value}),
      wilcox = apply(counts, 1, function(x){wilcox.test(x ~ condition)$p.value}),
      stringsAsFactors = FALSE
    ) %>% 
      mutate(FC = copd / health, 
             log2FC = log2(copd / health),
             fdrt = p.adjust(t, method = "fdr"),
             fdrwilcox = p.adjust(wilcox, method = "fdr"))
    
    save(counts, file = paste0(i,"_counts.Rdata"))
    save(DEG,    file = paste0(i,"_DEGs.Rdata"))
    write.csv(DEG, paste0(i,"_DEGs.csv"))
  }
  message(Sa(paste0("--> End ", i," DEA <--")))
}
