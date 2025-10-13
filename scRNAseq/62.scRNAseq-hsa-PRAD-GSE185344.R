require(Seurat)
require(dplyr)
require(ggplot2)
# load data---------------------------------------------------------------------
scobj.list <- readRDS("rawdat/GSE185344(rds)/GSE185344_PH_scRNA.final.rds")
annot <- data.table::fread("rawdat/GSE185344(rds)/GSE185344_PH_scRNA.rename_cluster.csv.gz")

names(scobj.list)
# [1] "rawobj"        "filter"        "markers"       "seurat_colors" "obj"  

class(scobj.list$markers)
class(scobj.list$rawobj)
class(scobj.list$seurat_colors)
class(scobj.list$obj)
names(scobj.list$filter)
# [1] "nFeature_cutoff_min" "nFeature_cutoff_max" "mt_cutoff"          
# [4] "nCount_cutoff"       "nCount_sd_cutoff"  

dim(scobj.list$rawobj)
# [1] 38224 62995
dim(scobj.list$obj)
# [1] 27107 57697

scobj <- scobj.list$obj
meta.data <- scobj@meta.data
meta.data$V1 <- rownames(meta.data)
meta.data <- meta.data %>% left_join(annot %>% dplyr::select(-seurat_clusters), by = "V1")
meta.data <- as.data.frame(meta.data)
rownames(meta.data) <- meta.data$V1
meta.data -> scobj@meta.data
# Visualization-----------------------------------------------------------------
DimPlot(scobj, reduction = "umap", 
        group.by = "renamed_cellactivity_clusters", 
        label = TRUE)

p_out <- lapply(c("AURKA","CDC20","KIF20A","PRR11"), 
                function(x){FeaturePlot(scobj, features = x) + 
                    scale_colour_gradientn(
                      colours = colorRampPalette(
                        c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                          '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                          '#6F003D','#56033F'))(1000))}) %>% 
  patchwork::wrap_plots(ncol = 2, nrow = 2)
# Visualization Manually--------------------------------------------------------
dat.umap <- scobj@reductions$umap@cell.embeddings %>% as.data.frame
dat.umap$V1 <- rownames(dat.umap)
dat.umap <- dat.umap %>% left_join(annot %>% dplyr::select(-seurat_clusters), by = "V1")

all(dat.umap$V1 == colnames(scobj@assays$RNA$data))
# [1] TRUE
dat.umap$AURKA <- scobj@assays$RNA$data['SYMBOL',] %>% unlist

p1 <- ggplot(dat.umap, aes(x = UMAP_1, y = UMAP_2, color = renamed_cellactivity_clusters)) +
  geom_point() + 
  labs(color = "") +
  theme_classic() +
  theme(panel.grid = element_blank(), legend)
p2 <- ggplot(dat.umap, aes(x = UMAP_1, y = UMAP_2, color = SYMBOL)) +
  geom_point(size = 0.1) + 
  scale_colour_gradientn(
    colours = colorRampPalette(
      c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
        '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
        '#6F003D','#56033F'))(1000)) +
  labs(color = "") +
  theme_classic() +
  theme(panel.grid = element_blank())

VlnPlot(scobj, features = c("SYMBOL","SYMBOL","SYMBOL","SYMBOL"), 
        assay = "RNA", layer = "data",
        group.by = "renamed_cellactivity_clusters", ncol = 4)

