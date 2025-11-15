rm(list = ls());gc()

library(dplyr)
library(stringr)
library(ggplot2)
library(biomaRt)
library(clusterProfiler)
# load data --------------------------------------------------------------------
count <- data.table::fread("GSE132143_HomoSapiens_Counts.txt.gz") %>% as.data.frame
colnames(count) <- colnames(count) %>% str_remove("\\.fastq\\.gz$")
samples <- colnames(count)[-1]

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

conversion_table <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "gene_biotype",
    # "ensembl_gene_id",
    "hgnc_symbol"),      
  filters = "ensembl_gene_id",
  values = count$V1,
  mart = human_mart
)

colnames(count)[1] <- "ensembl_gene_id"
count <- count %>% 
  left_join(conversion_table, by = "ensembl_gene_id") %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(hgnc_symbol != "") %>% 
  filter(complete.cases(hgnc_symbol))

symbol.dup <- count$hgnc_symbol[duplicated(count$hgnc_symbol)] %>% unique
count1 <- count %>% 
  filter(!(hgnc_symbol %in% symbol.dup)) %>% 
  dplyr::select(-gene_biotype, -ensembl_gene_id) %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("hgnc_symbol")
count2 <- count %>% 
  filter(hgnc_symbol %in% symbol.dup) %>% 
  tidyr::pivot_longer(samples) %>% 
  group_by(hgnc_symbol, name) %>% 
  mutate(sum = sum(value)) %>% 
  dplyr::select(-value) %>% 
  distinct(hgnc_symbol, name, .keep_all = TRUE) %>% 
  tidyr::pivot_wider(names_from = "name", values_from = "sum") %>% 
  filter(!duplicated(hgnc_symbol)) %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("hgnc_symbol")

count1  <- count1[,samples]
count2  <- count2[,samples]
count   <- rbind(count1, count2)

meta <- readxl::read_excel("GSE132143-GPL18573_series_matrix.xlsx") %>% arrange(group)
# DEA --------------------------------------------------------------------------
require(DESeq2)

count <- count[,meta$ID]
condition = factor(meta$group, levels = c("Healthy","MI","DCM"))
coldata <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~condition)
dds$condition <- relevel(dds$condition, ref = "Healthy") # 指定哪一组作为对照组
dds <- DESeq(dds) 
count.norm <- counts(dds, normalized=T) 
DEG.MI  <- results(dds, name="condition_MI_vs_Healthy", independentFiltering = FALSE) %>% as.data.frame
DEG.MI  <- na.omit(DEG.MI)
DEG.DCM <- results(dds, name="condition_DCM_vs_Healthy", independentFiltering = FALSE) %>% as.data.frame
DEG.DCM <- na.omit(DEG.DCM)

DEG.MI[c("IL1B","NEU1","TLR4","SIGLEC1"),]
#           baseMean log2FoldChange      lfcSE      stat     pvalue       padj
# IL1B      8.007674      0.4664615 0.53757895  0.867708 0.38555423 0.50877120
# NEU1    233.794529     -0.2194025 0.09097777 -2.411605 0.01588246 0.03195037
# TLR4    533.975034     -0.5719661 0.24090837 -2.374206 0.01758674 0.03492612
# SIGLEC1 244.532972      0.5992173 0.30200449  1.984134 0.04724093 0.08306530
# GSVA Glycosylation -----------------------------------------------------------
require(GSVA)

meta <- meta %>% filter(group != "DCM")
count.norm <- count.norm[,meta$ID]
dim(count.norm)
# [1] 19312    42

hallmarks <-
  c("GOBP_GLYCOSYLATION.v2025.1.Hs.gmt",                 
    "GOBP_PROTEIN_DEGLYCOSYLATION.v2025.1.Hs.gmt",       
    "GOBP_PROTEIN_N_LINKED_GLYCOSYLATION.v2025.1.Hs.gmt",
    "GOBP_PROTEIN_O_LINKED_GLYCOSYLATION.v2025.1.Hs.gmt")
res <- list()
n <- 1
for (i in hallmarks) {
  hallmark <- qusage::read.gmt(i) #返回的是表格
  gsvaPar <- ssgseaParam(exprData = as.matrix(count.norm), 
                         geneSets = hallmark,
                         normalize = TRUE)
  gsva_data <- gsva(gsvaPar, verbose = FALSE)
  res[[n]] <- gsva_data
  n <- n + 1
}

res.df <- rbind(res[[1]], res[[2]]) %>% 
  rbind(res[[3]]) %>% 
  rbind(res[[4]])
# visualization-----------------------------------------------------------------
annotation_rows <- data.frame(
  group = factor(meta$group, levels = c("Healthy","MI","DCM")),
  row.names = colnames(gsva_data)
)
pheatmap::pheatmap(res.df %>% t, 
                   clustering_method = "ward.D",
                   annotation_row = annotation_rows,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   scale = "column")

apply(res.df, 1, function(x){t.test(x ~ meta$group)})






