# data type: 基于质谱的蛋白组数据
# company: 青莲百奥
# sample: tissue from mice
# software: 


rm(list = ls());gc()
# load pkgs---------------------------------------------------------------------
require(dplyr)
require(stringr)
require(DESeq2)
require(clusterProfiler)
require(enrichplot)
# load data---------------------------------------------------------------------
# the data is intensity obtained by MS without normalization
# a series of methods for normalization can be found in Multi-Omics-Routine
# /Biostatistic/6. Normalization.R
dat <- readxl::read_excel("intensity.xlsx")
annot <- readxl::read_excel("pro2anno.xlsx")

boxplot(dat[,-1], outline = F, col = "red")

dim(dat)
# [1] 7910    7
# QC----------------------------------------------------------------------------
# Here, we calculate the sum of the proteins, who have the same name
dat <- dat %>% 
  filter(complete.cases(.)) %>% 
  left_join(annot %>% dplyr::select(Accession, `Gene Name`), by = "Accession")

dat.uni <- dat %>% 
  tidyr::pivot_longer(cols = c('U1','U2','U3','UF1','UF2','UF3'), 
                      names_to = "ID", values_to = "intensity") %>% 
  group_by(`Gene Name`, ID) %>% 
  mutate(sum = sum(intensity)) %>% 
  dplyr::select(-intensity) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = "ID", values_from = "sum") %>% 
  filter(!duplicated(`Gene Name`)) %>% 
  filter(complete.cases(`Gene Name`)) %>% 
  tibble::column_to_rownames("Gene Name")

dim(dat.uni)
# [1] 6032    7

dat.uni <- dat.uni[,-1]
# DESEQ2------------------------------------------------------------------------
data.cor <- as.data.frame(dat.uni)
group_list <- colnames(data.cor) %>% stringr::str_remove_all("[1-3]{1,2}$")
count <- data.cor
condition = factor(group_list, levels = c("U","UF"))
coldata <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(countData = round(count),
                              colData = coldata,
                              design = ~condition)
dds$condition<- relevel(dds$condition, ref = "U") # 指定哪一组作为对照组

# https://cloud.tencent.com/developer/article/2327198?from=15425
# results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
dds <- DESeq(dds)
DEG2_1 <- results(dds, name="condition_UF_vs_U", independentFiltering = FALSE) %>% as.data.frame
DEG2_1 <- na.omit(DEG2_1)
