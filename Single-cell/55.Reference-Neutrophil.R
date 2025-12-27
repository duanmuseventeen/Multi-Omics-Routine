# Single-cell transcriptome profiling reveals neutrophil heterogeneity in homeostasis and infection
# 2020
# Nat Immunol
# PMID: 32719519
# GSE137539

library(Seurat)
library(dplyr)
library(stringr)

# Load data---------------------------------------------------------------------
NKD <- data.table::fread("Xie/GSM4081555_NKD181202670.txt.gz") %>% 
  as.data.frame %>% tibble::column_to_rownames("Cell_Index") %>% t %>% 
  as.data.frame
RBD <- data.table::fread("Xie/GSM4473444_RBD00006.txt.gz") %>% 
  as.data.frame %>% tibble::column_to_rownames("V1")
HDY <- data.table::fread("Xie/GSM4473445_HDY.txt.gz") %>% 
  as.data.frame %>% tibble::column_to_rownames("V1")

lapply(list(NKD, RBD, HDY), dim)
lapply(list(NKD[,1:6], RBD[,1:6], HDY[,1:6]), head)
# filter gene symbols-----------------------------------------------------------
gene_symbol <- intersect(rownames(NKD), row.names(RBD)) %>% intersect(rownames(HDY))

NKD_comm <- NKD[gene_symbol,]
RBD_comm <- RBD[gene_symbol,]
HDY_comm <- HDY[gene_symbol,]

lapply(list(NKD_comm, RBD_comm, HDY_comm), dim)
# filter columns----------------------------------------------------------------
meta <- data.table::fread("Xie/donor_meta.txt")
all(meta$V1[meta$orig.ident == "D1"] %>% str_remove("^D[1-3]_") %in% colnames(NKD_comm))
all(meta$V1[meta$orig.ident == "D2"] %>% str_remove("^D[1-3]_") %in% colnames(RBD_comm))
all(meta$V1[meta$orig.ident == "D3"] %>% str_remove("^D[1-3]_") %in% colnames(HDY_comm))

all(meta$V1[meta$orig.ident == "D1"] %>% str_remove("^D[1-3]_") %in% colnames(RBD_comm))
all(meta$V1[meta$orig.ident == "D1"] %>% str_remove("^D[1-3]_") %in% colnames(HDY_comm))

all(meta$V1[meta$orig.ident == "D2"] %>% str_remove("^D[1-3]_") %in% colnames(NKD_comm))
all(meta$V1[meta$orig.ident == "D2"] %>% str_remove("^D[1-3]_") %in% colnames(HDY_comm))

all(meta$V1[meta$orig.ident == "D3"] %>% str_remove("^D[1-3]_") %in% colnames(NKD_comm))
all(meta$V1[meta$orig.ident == "D3"] %>% str_remove("^D[1-3]_") %in% colnames(RBD_comm))

# 未找到作者提供的对应关系，通过barcodes和样本细胞数推断:
# D1 - NKD
# D2 - RBD
# D3 - HDY 

D1 <- NKD_comm
colnames(D1) <- paste0("D1_", colnames(D1))
D1 <- D1[,meta$V1[meta$orig.ident == "D1"]]

D2 <- RBD_comm
colnames(D2) <- paste0("D2_", colnames(D2))
D2 <- D2[,meta$V1[meta$orig.ident == "D2"]]

D3 <- HDY_comm
colnames(D3) <- paste0("D3_", colnames(D3))
D3 <- D3[,meta$V1[meta$orig.ident == "D3"]]

all(c(colnames(D1), colnames(D2), colnames(D3)) == meta$V1)
# [1] TRUE
# Create Seurat Object----------------------------------------------------------
D1$SYMBOL <- rownames(D1)
D2$SYMBOL <- rownames(D2)
D3$SYMBOL <- rownames(D3)

raw.data <- D1 %>% 
  left_join(D2, by = "SYMBOL") %>% 
  left_join(D3, by = "SYMBOL") %>% 
  dplyr::select(-SYMBOL) %>% 
  as.data.frame
# dim(raw.data)
# [1] 19154 22030
rownames(raw.data) <- rownames(D1)
raw.data <- raw.data %>% dplyr::select(all_of(meta$V1))
# dim(raw.data)
# [1] 19154 22030
raw.data <- raw.data %>% 
  as.matrix %>% 
  Matrix::Matrix(sparse = TRUE)
scobj <- CreateSeuratObject(counts = raw.data, data = raw.data)
scobj <- AddMetaData(object = scobj, metadata = meta)

VlnPlot(scobj, features = 'nFeature_RNA')
VlnPlot(scobj, features = 'nCount_RNA')
VlnPlot(scobj, features = 'percent.mt')

scobj <- subset(scobj, subset = cell_type == "Neu")
table(scobj$cluster)
# Cont  hG5a  hG5b  hG5c 
# 31  4547  5738 10403

scobj <- subset(scobj, subset = cluster != "Cont")

Xie.2020 <- scobj 
save(Xie.2020, file = "C:/D/R project/Multi-Omics-Routine/scRNAseq/Neutrophil-Xie.2020.Rdata")


# Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species# 2019
# Immunity
# PMID: 30979687
# GSE127465
# Load data---------------------------------------------------------------------
dat <- Matrix::readMM("Zilionis/GSE127465_human_counts_normalized_54773x41861.mtx.gz")

class(data)
# [1] "function"
dim(dat)
# [1] 54773 41861
class(dat)
# [1] "dgTMatrix"
# attr(,"package")
# [1] "Matrix"

# add gene names
var <- data.table::fread('Zilionis/GSE127465_gene_names_human_41861.tsv.gz', header = FALSE)

# add per-cell metadata
meta<- data.table::fread('Zilionis/GSE127465_human_cell_metadata_54773x25.tsv.gz')
# Create Seurat Object----------------------------------------------------------
colnames(dat) <- var$V1
rownames(dat) <- paste0(meta$Library, "_", meta$Barcode)
scobj <- CreateSeuratObject(counts = t(dat), data = t(dat))
scobj <- AddMetaData(object = scobj, metadata = meta)
scobj$nCount_RNA <- scobj@meta.data$`Total counts`
scobj$perent.mt <- scobj@meta.data$`Percent counts from mitochondrial genes`
scobj$major_cell_type <- scobj@meta.data$`Major cell type`
scobj$minor_cell_type <- scobj@meta.data$`Minor subset`

VlnPlot(scobj, features = 'nFeature_RNA')
VlnPlot(scobj, features = 'nCount_RNA')
VlnPlot(scobj, features = 'perent.mt')

scobj <- subset(scobj, subset = major_cell_type %in% c("bNeutrophils", "tNeutrophils"))

table(scobj$major_cell_type)
# bNeutrophils tNeutrophils 
# 8911         2728 

table(scobj$minor_cell_type)
# bN1  bN2  bN3  bN4  bN5  bN6  tN1  tN2  tN3  tN4  tN5 
# 1331  646 2021 3387 1030  496  526   86  425  789  902

Zilions.2019 <- scobj 
save(Zilions.2019, file = "C:/D/R project/Multi-Omics-Routine/scRNAseq/Neutrophil-Zilions.2019.Rdata")
