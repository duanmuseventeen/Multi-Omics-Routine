# Transcriptional mediators of treatment resistance in lethal prostate cancer
# Nat Med
# 2021
# 33664492
# url: https://singlecell.broadinstitute.org/single_cell/study/SCP1244/transcriptional-mediators-of-treatment-resistance-in-lethal-prostate-cancer?genes=PRR11#study-summary

require(Seurat)
require(dplyr)
require(ggplot2)
require(Matrix)
require(patchwork)

# load data---------------------------------------------------------------------
cluster <- data.table::fread("scp_clustering.tsv")
metadata<- data.table::fread("scp_metadata.tsv")
nk_t_cluster <- data.table::fread("scp_nk_t_clustering.tsv")
tpm <- data.table::fread("scp_tpm.tsv.gz")
# V1    V2    V3    V4    V5    V6
# <char> <num> <num> <num> <num> <num>
# 1:       GENE     0     1     2     3     4
# 2: AL356585.2     0     0     0     0     0
# 3: CU638689.1     0     0     0     0     0
# 4: CU638689.2     0     0     0     0     0
# 5: CU634019.1     0     0     0     0     0
# 6: CU633906.1     0     0     0     0     0

# Create seurat object----------------------------------------------------------
tpm <- as.data.frame(tpm)
colnames(tpm) <- tpm[1,] %>% unlist
tpm <- tpm[-1,]
rownames(tpm) <- tpm$GENE
tpm <- tpm[,-1]

metadata <- metadata[-1,]
cluster <- cluster[-1,]
metadata <- metadata %>% 
  left_join(cluster, by = "NAME")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$NAME

tpm <- tpm %>% 
  as.matrix %>% 
  Matrix(sparse = TRUE)
scobj <- CreateSeuratObject(counts = tpm, data = tpm, meta.data = metadata, project = "Smart-seq2")

# Visualization-----------------------------------------------------------------
VlnPlot(scobj, features = c("AURKA","CDC20","KIF20A","PRR11"), 
        group.by = "cluster dominant cell type", ncol = 4)

dat.umap <- scobj@meta.data

all(rownames(dat.umap) == colnames(scobj@assays$RNA$data))
# [1] TRUE
dat.umap$SYMBOL <- scobj@assays$RNA$data['SYMBOL',] %>% unlist
dat.umap$X <- dat.umap$X %>% as.numeric
dat.umap$Y <- dat.umap$Y %>% as.numeric

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

p1 <- ggplot(dat.umap, aes(x = X, y = Y, color =`cluster dominant cell type`)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = colors_list) +
  labs(color = "") +
  theme_classic() +
  theme(panel.grid = element_blank())
    theme(panel.grid = element_blank())
p_feature <- ggplot(dat.umap, aes(x = X, y = Y, color = SYMBOL)) +
  geom_point(size = 0.1) + 
  scale_colour_gradientn(
    colours = colorRampPalette(
      c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
        '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
        '#6F003D','#56033F'))(1000)) +
  labs(color = "") +
  theme_classic() +
  theme(panel.grid = element_blank())

ggsave(p1, filename = "UMAP(33664492).pdf", width=6, height=5, units="in")
ggsave(p_feature, filename = "FeaturePlot(33664492).pdf", width=6, height=5, units="in")


