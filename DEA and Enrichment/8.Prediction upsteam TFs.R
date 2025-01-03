# screen all TFs
# because big size, rdata cannot be uploaded, please make in house db by yourself

require(dplyr)

search <- ""

# Homo species------------------------------------------------------------------
## 1. GTRD----------------------------------------------------------------------
# TF key: SYMBOL 
load("GTRD_hsa.db.Rdata")
## 2. KnockTF-------------------------------------------------------------------
# https://bio.liclab.net/KnockTFv2/download.php
# Detailed Info of KnockTF issues from the author
if(TRUE){
  C1: the sample ID of every T(co)F knockdown / knockout dataset
  
  C2: the TF that disrupted by different knockdown / knockout techniques
  
  C3: the target genes of Adnp
  
  C4: the average expression value of the disease sample
  
  C5: the average expression value of the normal sample
  
  C6: Fold Change, the multiples of mean expression values between disease and normal samples
  
  C7: the log 2 value of Fold Change
  
  C8: the sequencing results of differential expression of genes
  
  C9: the statistical significance calculated by limma
  
  C10: the FDR corrected P-value using the Benjamini-Hochberg method
  
  C11: The target gene is up-regulated or down-regulated
  
  C12: official symbol name of the target gene
  
  C13: Entrez ID of the target gene
  
  C14: Ensembl ID of the target gene
  
  C15: How many transcription factors are bound to the typical enhancer of the target gene
  
  C16: How many transcription factors are bound to the promoter region of the target gene
  
  C17: How many transcription factors are bound to the super-enhancer of the target gene
  
  C18: The name of a transcription factor that binds to the typical enhancer region of the target gene
  
  C19: The name of a transcription factor that binds to the promoter region of the target gene
  
  C20: The name of a transcription factor that binds to the super-enhancer region of the target gene
}

knockTF <- data.table::fread("differential expression of genes in all datasets(Human).txt")
## 3. hTFtarget-----------------------------------------------------------------
# http://bioinfo.life.hust.edu.cn/hTFtarget#!/
htTF <- data.table::fread("TF-Target-information.txt")
## 4. ENCODE--------------------------------------------------------------------
# https://maayanlab.cloud/Harmonizome/download
ENCODE <- data.table::fread("gene_attribute_edges.txt",skip = 1)
## 5. UCSC----------------------------------------------------------------------
## 6. Cistrome Db---------------------------------------------------------------
# http://cistrome.org/db/#/
# Application Failure
## 7. Iregulaton----------------------------------------------------------------
# please load result from Cytoscape::Iregulaton
## 8. JASPAR--------------------------------------------------------------------
# https://maayanlab.cloud/Harmonizome/download
JASPAR <- data.table::fread("gene_attribute_edges.txt",skip = 1)
## 9. HOCOMOCO------------------------------------------------------------------
# https://hocomoco11.autosome.org/
## 10.MotifMap------------------------------------------------------------------
# http://motifmap.ics.uci.edu/
# Download is not supported, if necessary, the data can be collected manually
## 11.TTRUST--------------------------------------------------------------------
# https://www.grnpedia.org/trrust/
# limited TFs in this db
TRUST <- data.table::fread("trrust_rawdata.human.tsv")
## 12.TF prediction db----------------------------------------------------------
# https://transcriptionfactor.org/?Home
## 13.Factorbook----------------------------------------------------------------
# http://v1.factorbook.org/mediawiki/index.php/Welcome_to_factorbook
## 14.Human Protein-DNA Interactome(hPDI)--------------------------------------- 
# http://bioinfo.wilmer.jhu.edu/PDI/
# cannot use search gene to find TF regulating it
## 15.TRANSFAC Matrix models http://gene-regulation.com/pub/databases.html------
# NOT FREE

# To find TF regulating search gene---------------------------------------------
res_gtrd <- gtrd_hsa %>% filter(`Gene symbol` == search) %>% dplyr::select(SYMBOL)
res_knockTF <- knockTF %>% filter(`Gene` == search) %>% dplyr::select(TF)
res_htTF <- htTF %>% filter(`target` == search) %>% dplyr::select(TF)
res_encode <- ENCODE %>% filter(`GeneSym` == search) %>% dplyr::select(TF)
res_jaspar <- JASPAR %>% filter(`GeneSym` == search) %>% dplyr::select(TF)
res_trust <- TRUST %>% filter(SYMBOL == search) %>% dplyr::select(TF)

# instersection
res_gtrd$SYMBOL %>% 
  intersect(res_knockTF$TF) %>%
  intersect(res_htTF$TF) %>%
  intersect(res_encode$TF) %>%
  intersect(res_jaspar$TF) %>%
  intersect(res_trust$TF)

# Mus Musculus------------------------------------------------------------------
## GTRD-------------------------------------------------------------------------
load("GTRD_mmu.db.Rdata")
## KnockTF----------------------------------------------------------------------
knockTF <- data.table::fread("differential expression of genes in all datasets(Mouse).txt")
## TRUST------------------------------------------------------------------------
TRUST <- data.table::fread("trrust_rawdata.mouse.tsv")

