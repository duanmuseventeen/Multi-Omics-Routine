rm(list = ls());gc()

setwd("")
# Load pkgs --------------------------------------------------------------------
library(limma)
library(stringr)
library(GEOquery)
library(enrichplot)
library(affy)
source("toolbox.R")
# load data --------------------------------------------------------------------
focus<- paste0("GAPDH")
data <- readxl::read_excel("GSE7392/GSE7392_expr.xlsx")
meta <- readxl::read_excel("GSE7392/GSE7392_meta.xlsx")
# annotate probe ---------------------------------------------------------------
GPL570 <- data.table::fread("GPL570-55999.txt", skip = 16)

samples   <- colnames(data)[-1]
data$mean <- rowMeans(data[,-1])
data.merge<- data %>% 
  left_join(GPL570, by = "ID") %>% 
  filter(complete.cases(`Gene Symbol`)) %>% 
  filter(`Gene Symbol` != "") %>% 
  filter(!str_detect(`Gene Symbol`, "///")) %>% 
  arrange(desc(mean)) %>% 
  filter(!duplicated(`Gene Symbol`))

dim(data.merge) # [1] 21655    47

focus %in% data.merge$`Gene Symbol`
# QC  --------------------------------------------------------------------------
expr <- data.merge %>% tibble::column_to_rownames("Gene Symbol") %>% 
  dplyr::select(all_of(samples))

meta <- meta %>% filter(geo_accession %in% colnames(expr))
expr <- expr[,meta$geo_accession]

range(expr)

# log2 transform
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  expr[which(expr <= 0)] <- NaN
  expr <- log2(expr + 1) 
  }

range(expr) # [1]  0.1219615 16.0838221

boxplot(expr)
# DEA --------------------------------------------------------------------------
group  <- meta$group %>% factor(levels = c(0, 1), labels = c("Normal","Fibrotic"))
design <- model.matrix(~ 0 + group)
rownames(design) <- colnames(expr)
colnames(design) <- levels(group)
cont.matrix <- makeContrasts(Fibrotic_vs_Normal = c('Fibrotic-Normal'), levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
degs <-topTable(fit2, coef = 'Fibrotic-Normal', number = Inf,  adjust.method = "BH", sort.by = "P")

degs[focus, ]
# Visualization ----------------------------------------------------------------
