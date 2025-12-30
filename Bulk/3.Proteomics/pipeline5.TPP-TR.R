setwd("TPP-TR/")

# Indolyl-3-propionic acid vs Control
# Indole-3-carboxylic acid vs Control
# Indole-3-lactic Acid     vs Control


rm(list = ls())
# load pkgs-------------------------------------------------------------------
require(dplyr)
require(stringr)
require(ggplot2)
require(TPP)

dat           <- readxl::read_excel("TPP.xlsx")
# preparation-------------------------------------------------------------------
# edit name---------------------------------------------------------------------
# omit

# QC & group--------------------------------------------------------------------
dat   <- dat %>% 
  mutate(gene_name = `Gene ID`,
         qssm = 0,
         qupm = 0) %>% 
  relocate(gene_name, qssm, qupm) %>% 
  tibble::column_to_rownames("Protein Name") %>% 
  dplyr::select(-`Protein Accession`, -`Gene ID`, -`Protein Description`, -`New Group avg`)
dat1 <- dat %>% dplyr::select(gene_name, qssm, qupm, contains("Control"))
dat2 <- dat %>% dplyr::select(gene_name, qssm, qupm, contains("ICA"))
dat3 <- dat %>% dplyr::select(gene_name, qssm, qupm, contains("IPA"))
dat4 <- dat %>% dplyr::select(gene_name, qssm, qupm, contains("ILA"))

colnames(dat1) <- colnames(dat1) %>% str_remove("Control")
colnames(dat2) <- colnames(dat2) %>% str_remove("ICA")
colnames(dat3) <- colnames(dat3) %>% str_remove("IPA")
colnames(dat4) <- colnames(dat4) %>% str_remove("ILA")

dat1           <- dat1 %>% dplyr::select(gene_name, qssm, qupm, `67`,`63`,`56`,`50`,`44`,`37`)
dat2           <- dat2 %>% dplyr::select(gene_name, qssm, qupm, `67`,`63`,`56`,`50`,`44`,`37`)
dat3           <- dat3 %>% dplyr::select(gene_name, qssm, qupm, `67`,`63`,`56`,`50`,`44`,`37`)
dat4           <- dat4 %>% dplyr::select(gene_name, qssm, qupm, `67`,`63`,`56`,`50`,`44`,`37`)

myfc <- function(data){
  data      <- data %>% filter(`37` != 0)
  ref       <- data$`37`
  data$`37` <- data$`37`/ ref
  data$`44` <- data$`44`/ ref
  data$`50` <- data$`50`/ ref
  data$`56` <- data$`56`/ ref
  data$`63` <- data$`63`/ ref
  data$`67` <- data$`67`/ ref
  colnames(data)[4:ncol(data)] <- paste0("rel_fc_", colnames(data)[4:ncol(data)])
  return(data)
}

dat1 <- myfc(data = dat1)
dat2 <- myfc(data = dat2) # ICA
dat3 <- myfc(data = dat3) # IPA
dat4 <- myfc(data = dat4) # ILA
# 2.ICA---------------------------------------------------------------------------
commonname <- intersect(rownames(dat1), rownames(dat2))
dat1.ica   <- dat1[commonname,]
dat2.ica   <- dat2[commonname,]

# 首先定义一个新的归一化标准，适应你的 6 个温度点
myNormReqs <- tpptrDefaultNormReqs()
myNormReqs$nFCHoldout <- 6  # 将要求的有效 Fold Change 数量改为 6

TRresults <- analyzeTPPTR(
  configTable = data.frame(
    Experiment = c("Control","ICA"),
    Condition  = c("Vehicle","Treatment"),
    ComparisonVT1 = c("x", "x"),
    `37` = rep(37, 2),
    `44` = rep(44, 2),
    `50` = rep(50, 2),
    `56` = rep(56, 2),
    `63` = rep(63, 2),
    `67` = rep(67, 2),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  normalize = FALSE,
  methods = "meltcurvefit",
  data = list(Control = dat1.ica, ICA = dat2.ica),
  nCores = 8,
  resultPath = "2.ICA",
  normReqs = myNormReqs,  # 在这里传入你修改后的标准
  plotCurves = FALSE)
