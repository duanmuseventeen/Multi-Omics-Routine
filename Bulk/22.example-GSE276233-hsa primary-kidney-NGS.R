require(DESeq2)
require(stringr)
require(dplyr)

datList <- lapply(dir(), function(x){
  tmp <- data.table::fread(x) %>% 
    filter(BioType %in% c(#"lincRNA", "macro_lncRNA", "miRNA", 
      "protein_coding")) %>% 
    dplyr::select(`2ndStrandCount`, GeneName)
  colnames(tmp) <- c(x, "SYMBOL")
  return(tmp)
})

for (i in 1:12) {
  print(identical(datList[[1]]$SYMBOL, datList[[i]]$SYMBOL))
}

dat <- do.call(bind_cols, datList)
dat$`GENE SYMBOL` <- dat$`SYMBOL...2`
dat <- dat %>% dplyr::select(-starts_with("SYMBOL"))
colnames(dat) <- str_remove(colnames(dat), ".tab")
dat <- as.data.frame(dat)

all(dat$`GENE SYMBOL` == datList[[3]]$SYMBOL)

dat.uniq <- dat %>% 
  tidyr::pivot_longer(names_to = "sample", values_to = "expr", 
                      cols = colnames(dat)[1:12]) %>% 
  group_by(`GENE SYMBOL`, sample) %>% 
  mutate(sum = sum(expr)) %>% 
  filter(!duplicated(paste0(sample, `GENE SYMBOL`))) %>% 
  dplyr::select(-expr) %>% 
  tidyr::pivot_wider(names_from = "sample", values_from = "sum") %>% 
  as.data.frame

dim(dat)
# [1] 29222    13
dim(dat.uniq)
# [1] 29192    13
sum(duplicated(dat$`GENE SYMBOL`))
# [1] 30

dat.uniq <- dat.uniq %>% tibble::column_to_rownames("GENE SYMBOL")

meta <- readxl::read_excel("../GSE276233_series_matrix/GSE276233_series_matrix.xlsx")
all(meta$geo_accession == colnames(dat.uniq))
# [1] TRUE

count <- dat.uniq
condition = factor(meta$group, levels = c("PBS","CSRP3"))
coldata <- data.frame(row.names = colnames(count), condition)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~condition)
dds$condition <- relevel(dds$condition, ref = "PBS") # 指定哪一组作为对照组
dds <- DESeq(dds)  
DEG <- results(dds, name="condition_CSRP3_vs_PBS", independentFiltering = FALSE) %>% as.data.frame
DEG <- na.omit(DEG)


















