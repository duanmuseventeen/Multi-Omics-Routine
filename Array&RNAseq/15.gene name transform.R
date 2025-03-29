# https://www.jianshu.com/p/abf5138e3757

library(dplyr)
library(stringr)
library(clusterProfiler)

# Section I: genes within different species-------------------------------------
hpa <- data.table::fread("DataBase/ProteinAtlas/proteinatlas/proteinatlas.tsv")

# BiocManager::install("homologene")
library(homologene)
# mouse: 10090
# human: 9606
homologene(hpa$Gene, inTax = 9606, outTax = 10090)
# Section II: genes within different sources------------------------------------
hpa <- data.table::fread("DataBase/ProteinAtlas/proteinatlas/proteinatlas.tsv")

annot <- bitr(hpa$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")



