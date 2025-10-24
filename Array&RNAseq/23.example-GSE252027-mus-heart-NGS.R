rm(list = ls());gc()

library(dplyr)
library(stringr)
library(ggplot2)
library(biomaRt)

if(CM){
  # load data --------------------------------------------------------------------
  setwd("CM/")
  datList <- lapply(dir(), function(x){
    tmp <- data.table::fread(x) %>% 
      dplyr::select(Name, TPM, NumReads)
    colnames(tmp)[2:3] <- paste0(str_remove(x, "\\.quant\\.genes\\.sf\\.gz$"), "_",
                                 colnames(tmp)[2:3])
    tmp <- tmp %>% as.data.frame %>% tibble::column_to_rownames("Name")
    return(tmp)
  })
  
  for (i in 2:24) {
    print(all(rownames(datList[[1]]) == rownames(datList[[2]])))
  }
  
  dat <- do.call(bind_cols, datList)
  # annot ------------------------------------------------------------------------
  # mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  conversion_table <- getBM(
    attributes = c(
      "ensembl_transcript_id", # 保持输入 ID，用于匹配
      "mgi_symbol",                    # 对应的基因符号（Gene Symbol）
      "external_gene_name"             # 外部基因名（有时与mgi_symbol不同，作为备份）
    ),
    filters = "ensembl_transcript_id",
    values = str_remove(rownames(dat), "\\.\\d+$"),
    mart = mouse_mart
  )
  
  dat.annot <- dat %>% 
    tibble::rownames_to_column("ensembl_transcript_id") %>% 
    mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.\\d+$")) %>% 
    left_join(conversion_table, by = "ensembl_transcript_id")
  
  # filter----
  dat.annot <- dat.annot %>% 
    filter(mgi_symbol != "") %>% 
    filter(complete.cases(mgi_symbol))
  
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "Gm\\d+$")),]
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "[0-9A-Z]+Rik$")),]
  
  dat.count <- dat.annot %>% 
    dplyr::select(contains("_NumReads")) %>% 
    mutate(symbol = dat.annot$mgi_symbol)
  
  dat.count.uniq <- dat.count %>% 
    tidyr::pivot_longer(cols = colnames(dat.count)[1:24]) %>% 
    group_by(name, symbol) %>%
    mutate(sum = round(sum(value))) %>% 
    filter(!duplicated(paste0(name, symbol))) %>% 
    dplyr::select(-value) %>% 
    tidyr::pivot_wider(names_from = "name", values_from = "sum")
  
  dat.count.uniq <- dat.count.uniq %>% 
    as.data.frame %>% 
    tibble::column_to_rownames("symbol") %>% 
    as.data.frame()
  
  # dat.tmp <- dat.annot %>% 
  #   dplyr::select(contains("_TPM")) %>% 
  #   mutate(symbol = dat.annot$mgi_symbol)
  
  # DEA --------------------------------------------------------------------------
  require(DESeq2)
  
  count <- dat.count.uniq
  condition = factor(colnames(dat.count.uniq) %>% 
                       str_remove("_CM_\\d_NumReads$") %>% 
                       str_remove("GSM\\d+_"), 
                     levels = c("ctrl","d1","d3","d7","d14","d28"))
  coldata <- data.frame(row.names = colnames(count), condition)
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = coldata,
                                design = ~condition)
  dds$condition <- relevel(dds$condition, ref = "ctrl") # 指定哪一组作为对照组
  dds <- DESeq(dds)  
  DEG1 <- results(dds, name="condition_d1_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG3 <- results(dds, name="condition_d3_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG7 <- results(dds, name="condition_d7_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG14<- results(dds, name="condition_d14_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG28<- results(dds, name="condition_d28_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  
  DEG1 <- na.omit(DEG1)
  DEG3 <- na.omit(DEG3)
  DEG7 <- na.omit(DEG7)
  DEG14<- na.omit(DEG14)
  DEG28<- na.omit(DEG28)
}

if(EC){
  # load data --------------------------------------------------------------------
  setwd("EC/")
  datList <- lapply(dir(), function(x){
    tmp <- data.table::fread(x) %>% 
      dplyr::select(Name, TPM, NumReads)
    colnames(tmp)[2:3] <- paste0(str_remove(x, "\\.quant\\.genes\\.sf\\.gz$"), "_",
                                 colnames(tmp)[2:3])
    tmp <- tmp %>% as.data.frame %>% tibble::column_to_rownames("Name")
    return(tmp)
  })
  
  for (i in 2:24) {
    print(all(rownames(datList[[1]]) == rownames(datList[[2]])))
  }
  
  dat <- do.call(bind_cols, datList)
  # annot ------------------------------------------------------------------------
  mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  conversion_table <- getBM(
    attributes = c(
      "ensembl_transcript_id", # 保持输入 ID，用于匹配
      "mgi_symbol",                    # 对应的基因符号（Gene Symbol）
      "external_gene_name"             # 外部基因名（有时与mgi_symbol不同，作为备份）
    ),
    filters = "ensembl_transcript_id",
    values = str_remove(rownames(dat), "\\.\\d+$"),
    mart = mouse_mart
  )
  
  dat.annot <- dat %>% 
    tibble::rownames_to_column("ensembl_transcript_id") %>% 
    mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.\\d+$")) %>% 
    left_join(conversion_table, by = "ensembl_transcript_id")
  
  # filter----
  dat.annot <- dat.annot %>% 
    filter(mgi_symbol != "") %>% 
    filter(complete.cases(mgi_symbol))
  
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "Gm\\d+$")),]
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "[0-9A-Z]+Rik$")),]
  
  dat.count <- dat.annot %>% 
    dplyr::select(contains("_NumReads")) %>% 
    mutate(symbol = dat.annot$mgi_symbol)
  
  dat.count.uniq <- dat.count %>% 
    tidyr::pivot_longer(cols = colnames(dat.count)[1:24]) %>% 
    group_by(name, symbol) %>%
    mutate(sum = round(sum(value))) %>% 
    filter(!duplicated(paste0(name, symbol))) %>% 
    dplyr::select(-value) %>% 
    tidyr::pivot_wider(names_from = "name", values_from = "sum")
  
  dat.count.uniq <- dat.count.uniq %>% 
    as.data.frame %>% 
    tibble::column_to_rownames("symbol") %>% 
    as.data.frame()
  
  # dat.tmp <- dat.annot %>% 
  #   dplyr::select(contains("_TPM")) %>% 
  #   mutate(symbol = dat.annot$mgi_symbol)
  
  # DEA --------------------------------------------------------------------------
  require(DESeq2)
  
  count <- dat.count.uniq
  condition = factor(colnames(dat.count.uniq) %>% 
                       str_remove("_EC_\\d_NumReads$") %>% 
                       str_remove("GSM\\d+_"), 
                     levels = c("ctrl","d1","d3","d7","d14","d28"))
  coldata <- data.frame(row.names = colnames(count), condition)
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = coldata,
                                design = ~condition)
  dds$condition <- relevel(dds$condition, ref = "ctrl") # 指定哪一组作为对照组
  dds <- DESeq(dds)  
  DEG1 <- results(dds, name="condition_d1_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG3 <- results(dds, name="condition_d3_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG7 <- results(dds, name="condition_d7_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG14<- results(dds, name="condition_d14_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG28<- results(dds, name="condition_d28_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  
  DEG1 <- na.omit(DEG1)
  DEG3 <- na.omit(DEG3)
  DEG7 <- na.omit(DEG7)
  DEG14<- na.omit(DEG14)
  DEG28<- na.omit(DEG28)
}

if(FB){
  # load data --------------------------------------------------------------------
  setwd("FB/")
  datList <- lapply(dir(), function(x){
    tmp <- data.table::fread(x) %>% 
      dplyr::select(Name, TPM, NumReads)
    colnames(tmp)[2:3] <- paste0(str_remove(x, "\\.quant\\.genes\\.sf\\.gz$"), "_",
                                 colnames(tmp)[2:3])
    tmp <- tmp %>% as.data.frame %>% tibble::column_to_rownames("Name")
    return(tmp)
  })
  
  for (i in 2:24) {
    print(all(rownames(datList[[1]]) == rownames(datList[[2]])))
  }
  
  dat <- do.call(bind_cols, datList)
  # annot ------------------------------------------------------------------------
  mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  conversion_table <- getBM(
    attributes = c(
      "ensembl_transcript_id", # 保持输入 ID，用于匹配
      "mgi_symbol",                    # 对应的基因符号（Gene Symbol）
      "external_gene_name"             # 外部基因名（有时与mgi_symbol不同，作为备份）
    ),
    filters = "ensembl_transcript_id",
    values = str_remove(rownames(dat), "\\.\\d+$"),
    mart = mouse_mart
  )
  
  dat.annot <- dat %>% 
    tibble::rownames_to_column("ensembl_transcript_id") %>% 
    mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.\\d+$")) %>% 
    left_join(conversion_table, by = "ensembl_transcript_id")
  
  # filter----
  dat.annot <- dat.annot %>% 
    filter(mgi_symbol != "") %>% 
    filter(complete.cases(mgi_symbol))
  
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "Gm\\d+$")),]
  dat.annot <- dat.annot[is.na(str_extract(dat.annot$mgi_symbol, "[0-9A-Z]+Rik$")),]
  
  dat.count <- dat.annot %>% 
    dplyr::select(contains("_NumReads")) %>% 
    mutate(symbol = dat.annot$mgi_symbol)
  
  dat.count.uniq <- dat.count %>% 
    tidyr::pivot_longer(cols = colnames(dat.count)[1:24]) %>% 
    group_by(name, symbol) %>%
    mutate(sum = round(sum(value))) %>% 
    filter(!duplicated(paste0(name, symbol))) %>% 
    dplyr::select(-value) %>% 
    tidyr::pivot_wider(names_from = "name", values_from = "sum")
  
  dat.count.uniq <- dat.count.uniq %>% 
    as.data.frame %>% 
    tibble::column_to_rownames("symbol") %>% 
    as.data.frame()
  
  # dat.tmp <- dat.annot %>% 
  #   dplyr::select(contains("_TPM")) %>% 
  #   mutate(symbol = dat.annot$mgi_symbol)
  
  # DEA --------------------------------------------------------------------------
  require(DESeq2)
  
  count <- dat.count.uniq
  condition = factor(colnames(dat.count.uniq) %>% 
                       str_remove("_FB_\\d_NumReads$") %>% 
                       str_remove("GSM\\d+_"), 
                     levels = c("ctrl","d1","d3","d7","d14","d28"))
  coldata <- data.frame(row.names = colnames(count), condition)
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = coldata,
                                design = ~condition)
  dds$condition <- relevel(dds$condition, ref = "ctrl") # 指定哪一组作为对照组
  dds <- DESeq(dds)  
  DEG1 <- results(dds, name="condition_d1_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG3 <- results(dds, name="condition_d3_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG7 <- results(dds, name="condition_d7_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG14<- results(dds, name="condition_d14_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  DEG28<- results(dds, name="condition_d28_vs_ctrl", independentFiltering = FALSE) %>% as.data.frame
  
  DEG1 <- na.omit(DEG1)
  DEG3 <- na.omit(DEG3)
  DEG7 <- na.omit(DEG7)
  DEG14<- na.omit(DEG14)
  DEG28<- na.omit(DEG28)
}



