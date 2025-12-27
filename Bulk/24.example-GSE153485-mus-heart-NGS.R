rm(list = ls());gc()

library(dplyr)
library(stringr)
library(ggplot2)
# load data --------------------------------------------------------------------
dat <- data.table::fread("GSE153485_raw_counts.txt.gz")
meta <- data.table::fread("GSE153485_metadata_summary.txt.gz")

sum(duplicated(dat$`Gene name`))

dat <- dat[is.na(str_extract(dat$`Gene name`, "Gm\\d+$")),]
dat <- dat[is.na(str_extract(dat$`Gene name`, "[0-9A-Z]+Rik$")),]
dat <- dat %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("Gene name")
dat <- round(dat)

dat <- dat[!(apply(dat, 1, function(x){any(x) == 0})),]
# DEA --------------------------------------------------------------------------
require(DESeq2)

# heart----
meta.heart <- meta %>% filter(Tissue == "heart")
# heart 6 h----
meta.heart.6 <- meta.heart %>% filter(Time == "6h")

count <- dat[,meta.heart.6$`NGI ID`]
condition = factor(meta.heart.6$Treatment, 
                   levels = c("SHAM","MI"))
coldata <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~condition)
dds$condition <- relevel(dds$condition, ref = "SHAM") # 指定哪一组作为对照组
dds <- DESeq(dds)  
DEG.heart.6 <- results(dds, name="condition_MI_vs_SHAM", independentFiltering = FALSE) %>% as.data.frame
DEG.heart.6 <- na.omit(DEG.heart.6)

# heart 24 h----
meta.heart.24<- meta.heart %>% filter(Time == "24h")

count <- dat[,meta.heart.24$`NGI ID`]
condition = factor(meta.heart.24$Treatment, 
                   levels = c("SHAM","MI"))
coldata <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~condition)
dds$condition <- relevel(dds$condition, ref = "SHAM") # 指定哪一组作为对照组
dds <- DESeq(dds)  
DEG.heart.24 <- results(dds, name="condition_MI_vs_SHAM", independentFiltering = FALSE) %>% as.data.frame
DEG.heart.24 <- na.omit(DEG.heart.24)

# GSVA -------------------------------------------------------------------------
require(GSVA)

# 6 h----
meta.heart.6 <- meta.heart.6 %>% arrange(desc(Treatment))

gsva.6 <- function(pathway){
  hallmark <- qusage::read.gmt(pathway)
  gsvaPar <- ssgseaParam(exprData = as.matrix(dat[,meta.heart.6$`NGI ID`]), 
                         geneSets = hallmark,
                         normalize = TRUE)
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE) %>% t
  return(gsva_matrix)
}
gsva.matrixList <- lapply(list(
  GOBP_CHOLINE_METABOLIC_PROCESS = "GOBP_CHOLINE_METABOLIC_PROCESS.v2025.1.Mm.gmt",
  GOBP_CHOLINE_TRANSPORT = "GOBP_CHOLINE_TRANSPORT.v2025.1.Mm.gmt",
  GOMF_CHOLINE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY = "GOMF_CHOLINE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY.v2025.1.Mm.gmt",
  REACTOME_CHOLINE_CATABOLISM = "REACTOME_CHOLINE_CATABOLISM.v2025.1.Mm.gmt"
), gsva.6)

gsva.heart.6 <- do.call(bind_cols, gsva.matrixList) %>% as.data.frame
rownames(gsva.heart.6) <- gsva.matrixList$GOBP_CHOLINE_METABOLIC_PROCESS %>% rownames
for (i in 1:4) {
  print(all(rownames(gsva.heart.6) == rownames(gsva.matrixList[[i]])))
}

annotation_rows <- data.frame(
  Group = meta.heart.6$Treatment %>% factor(levels = c("SHAM","MI")),
  row.names = meta.heart.6$`NGI ID`
)
pheatmap::pheatmap(gsva.heart.6, 
                   clustering_method = "ward.D",
                   annotation_row = annotation_rows,
                   cluster_rows = F, scale = "column")
apply(gsva.heart.6, 2, function(x){
  t.test(x ~ meta.heart.6$Treatment)
})

# 24 h----
meta.heart.24 <- meta.heart.24 %>% arrange(desc(Treatment))

gsva.24 <- function(pathway){
  hallmark <- qusage::read.gmt(pathway)
  gsvaPar <- ssgseaParam(exprData = as.matrix(dat[,meta.heart.24$`NGI ID`]), 
                         geneSets = hallmark,
                         normalize = TRUE)
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE) %>% t
  return(gsva_matrix)
}
gsva.matrixList <- lapply(list(
  GOBP_CHOLINE_METABOLIC_PROCESS = "GOBP_CHOLINE_METABOLIC_PROCESS.v2025.1.Mm.gmt",
  GOBP_CHOLINE_TRANSPORT = "GOBP_CHOLINE_TRANSPORT.v2025.1.Mm.gmt",
  GOMF_CHOLINE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY = "GOMF_CHOLINE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY.v2025.1.Mm.gmt",
  REACTOME_CHOLINE_CATABOLISM = "REACTOME_CHOLINE_CATABOLISM.v2025.1.Mm.gmt"
), gsva.24)

gsva.heart.24 <- do.call(bind_cols, gsva.matrixList) %>% as.data.frame
rownames(gsva.heart.24) <- gsva.matrixList$GOBP_CHOLINE_METABOLIC_PROCESS %>% rownames
for (i in 1:4) {
  print(all(rownames(gsva.heart.24) == rownames(gsva.matrixList[[i]])))
}

annotation_rows <- data.frame(
  Group = meta.heart.24$Treatment %>% factor(levels = c("SHAM","MI")),
  row.names = meta.heart.24$`NGI ID`
)
pheatmap::pheatmap(gsva.heart.24, 
                   clustering_method = "ward.D",
                   annotation_row = annotation_rows,
                   cluster_rows = F, scale = "column")

apply(gsva.heart.24, 2, function(x){
  t.test(x ~ meta.heart.24$Treatment)
})

save(gsva.heart.6, gsva.heart.24, file = "GSE153485_GSVA.Rdata")












