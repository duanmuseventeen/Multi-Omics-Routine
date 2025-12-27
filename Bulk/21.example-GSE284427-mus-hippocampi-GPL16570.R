# GSE284427 --------------------------------------------------------------------
# NEU1 wild or deficient 
# GPL16570

# load pkgs --------------------------------------------------------------------
# library(affy)
library(oligo)
library(dplyr)
library(stringr)
library(limma)
# load data --------------------------------------------------------------------
celFiles <- list.celfiles('GSE284427_RAW/', full.names=TRUE)
rawData <- read.celfiles(celFiles)

meta <- readxl::read_excel("GSE284427_series_matrix/meta.xlsx")
# PreProccess --------------------------------------------------------------------
# MAplot(rawData, pairs=TRUE)
# image(rawData[1])
# pmSeq <- pmSequence(rawData)

boxplot(rawData, target = "core")
hist(rawData, target = "core")

ppData <- rma(rawData)
boxplot(ppData, target = "core")

expr <- exprs(ppData) %>% as.data.frame
colnames(expr) <- str_extract(colnames(expr), "^GSM[0-9]+")
# AnnoProbe --------------------------------------------------------------------
annot <- data.table::fread("GSE284427_series_matrix/GPL16570-1802.txt")

myannot <- annot %>% 
  dplyr::select(ID, gene_assignment) %>% 
  filter(gene_assignment != "---") %>% 
  filter(gene_assignment != "") %>% 
  mutate(symbol = NA)

for (i in 1:nrow(myannot)) {
  myannot$symbol[i] <- unlist(str_split(myannot$gene_assignment[i],"//"))[2]
  myannot$symbol[i] <- str_remove_all(myannot$symbol[i], " ")
  print(i)
}

expr$ID <- rownames(expr)
dat <- expr %>%
  filter(ID %in% myannot$ID) %>% 
  left_join(myannot %>% dplyr::select(-gene_assignment), by = "ID") %>% 
  dplyr::select(-ID) %>% 
  tidyr::pivot_longer(names_to = "sample", values_to = "expr", cols = colnames(expr)[1:12]) %>% 
  group_by(symbol, sample) %>% 
  mutate(mean = mean(expr)) %>% 
  filter(!duplicated(paste0(sample, symbol))) %>% 
  dplyr::select(-expr) %>% 
  tidyr::pivot_wider(names_from = "sample", values_from = "mean")

dim(dat)
# [1] 24647    13
dim(myannot)
# [1] 27037     3
sum(duplicated(myannot$symbol))
# [1] 2390

dat <- dat %>% tibble::column_to_rownames("symbol")
# limma --------------------------------------------------------------------
dat.1m <- dat %>% dplyr::select(all_of(meta$geo_accession[meta$`age (month)` == 1]))
dat.5m <- dat %>% dplyr::select(all_of(meta$geo_accession[meta$`age (month)` == 5]))

meta$genotype <- str_replace(meta$genotype, "-", "_")
meta$genotype <- str_replace(meta$genotype, " ", "_")

# 1 month----
meta.1m <- meta %>% filter(`age (month)` == 1)
group.1m = factor(meta.1m$genotype, levels = c("Wild_type","Neu1_deficient"))
design <- model.matrix(~0+group.1m)
colnames(design) <- c("Wild_type","Neu1_deficient")

fit <- lmFit(dat.1m, design)
cont.matrix <- makeContrasts(contrasts = 'Neu1_deficient-Wild_type', levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

limma.1m = topTable(fit2, coef = 'Neu1_deficient-Wild_type', n = Inf)
limma.1m = na.omit(limma.1m)

# 5 month----
meta.5m <- meta %>% filter(`age (month)` == 5)
group.5m = factor(meta.5m$genotype, levels = c("Wild_type","Neu1_deficient"))
design <- model.matrix(~0+group.5m)
colnames(design) <- c("Wild_type","Neu1_deficient")

fit <- lmFit(dat.5m, design)
cont.matrix <- makeContrasts(contrasts = 'Neu1_deficient-Wild_type', levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

limma.5m = topTable(fit2, coef = 'Neu1_deficient-Wild_type', n = Inf)
limma.5m = na.omit(limma.5m)







