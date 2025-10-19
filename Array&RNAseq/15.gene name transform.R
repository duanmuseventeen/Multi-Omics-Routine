# https://www.jianshu.com/p/abf5138e3757
# https://blog.csdn.net/weixin_42655515/article/details/111167079
# https://zhuanlan.zhihu.com/p/593447503

library(dplyr)
library(stringr)
library(clusterProfiler)
library(Seurat)

# Section I: genes within different species-------------------------------------
# Method 1
hpa <- data.table::fread("DataBase/ProteinAtlas/proteinatlas/proteinatlas.tsv")

# BiocManager::install("homologene")
library(homologene)
# mouse: 10090
# human: 9606
homologene(hpa$Gene, inTax = 9606, outTax = 10090)

# Method 2
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

library(biomaRt)
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

conversion_table <- getLDS(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = c(s.genes, g2m.genes),
  mart = human_mart,
  attributesL = c("mgi_symbol"),
  martL = mouse_mart,
  uniqueRows = TRUE 
)

# Section II: genes within different sources------------------------------------
hpa <- data.table::fread("DataBase/ProteinAtlas/proteinatlas/proteinatlas.tsv")

annot <- bitr(hpa$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")




