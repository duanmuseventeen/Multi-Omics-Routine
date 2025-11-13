# ref: https://rpubs.com/LiYumei/806213

rm(list = ls());gc()
# load pkgs --------------------------------------------------------------------
require(dplyr)
require(stringr)
require(ggplot2)
require(survminer)
require(survival)
require(ggpubr)
require(edgeR)
# load data --------------------------------------------------------------------
load("TCGA-PRAD.Rdata")

# Preparation ------------------------------------------------------------------
count <- PRAD_counts %>% 
  filter(gene_type %in% c("protein_coding")) %>% 
  dplyr::select(-gene_id, -gene_type) %>% 
  as.data.frame

samples <- colnames(PRAD_counts)[-c(1:3)]
symbol.dup <- count$gene_name[duplicated(count$gene_name)] %>% unique
count1 <- count %>% 
  filter(!(gene_name %in% symbol.dup)) %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("gene_name")
count2 <- count %>% 
  filter(gene_name %in% symbol.dup) %>% 
  tidyr::pivot_longer(samples) %>% 
  group_by(gene_name, name) %>% 
  mutate(sum = sum(value)) %>% 
  dplyr::select(-value) %>% 
  distinct(gene_name, name, .keep_all = TRUE) %>% 
  tidyr::pivot_wider(names_from = "name", values_from = "sum") %>% 
  filter(!duplicated(gene_name)) %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("gene_name")

count1  <- count1[,samples]
count2  <- count2[,samples]
count   <- rbind(count1, count2)

group   <- str_sub(colnames(count), 14, 15)
count   <- count[, group %in% c("01", "11")]
group   <- group[group %in% c("01", "11")]  
group   <- ifelse(group == "11", 0, 1)
# DEA --------------------------------------------------------------------------
y <- DGEList(counts=count, group=group)
##Remove rows consistently have zero or very low counts
# keep <- filterByExpr(y)
# y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
count.n <- calcNormFactors(y, method="TMM") %>% cpm %>% as.data.frame

res <- data.frame(
  SYMBOL = rownames(count.n),
  norm.mean = count.n[,group == 0] %>% rowMeans,
  norm.sd   = apply(count.n[,group == 0], 1, sd),
  tumo.mean = count.n[,group == 1] %>% rowMeans,
  tumo.sd   = apply(count.n[,group == 1], 1, sd),
  FC = NA,
  wilcox.p = apply(count.n, 1, function(x){wilcox.test(x ~ group)$p.value})) %>%
  mutate(FC = tumo.mean / norm.mean,
         fdr = p.adjust(wilcox.p, "fdr"),
         bon = p.adjust(wilcox.p, "bonferroni")) %>%
  mutate(log2FC = log2(FC))

# wilcox.test(unlist(count.n["ACTB",]) ~ group)
# res[res$SYMBOL == "ACTB",]

save(res, file = "TCGA-PRAD(cpm-wilcox).Rdata")
write.csv(res, "TCGA-PRAD(cpm-wilcox).csv")