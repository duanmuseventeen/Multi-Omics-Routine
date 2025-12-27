rm(list = ls());gc()

setwd("~/SCP2155/expression/")
set.seed(1011)
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(scDblFinder)
library(DoubletFinder)
library(Augur)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(patchwork)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(GSVA)
# Load data---------------------------------------------------------------------
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB
set.seed(1011)

COPD  -   chronic obstructive pulmonary disease
IPF   -   idiopathic pulmonary fibrosis 
SSc   -   Systemic Sclerosis

ATI/II-   type I/II pneumocyte
  
  
colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

scobj <- "SCP2155" %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "SCP2155",
    assay = "RNA")

meta.data = data.table::fread("~/SCP2155/2022_SI_human_lung_allcells_metadata_file_updated.txt")
meta.data <- as.data.frame(meta.data)
rownames(meta.data) <- meta.data$NAME

all(rownames(scobj@meta.data) == rownames(meta.data))

scobj@meta.data <- meta.data
save(scobj, file = "SCP2155_processed.Rdata")

cluster <- data.table::fread("~/SCP2155/2022_SI_human_lung_allcells_cluster_file.txt") %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("NAME") %>% 
  as.matrix
colnames(cluster) <- c("UMAP_1", "UMAP_2")
umap <- CreateDimReducObject(
  embeddings = cluster, 
  key = "UMAP_", # 降维名称的前缀
  assay = DefaultAssay(scobj)
)
scobj@reductions$umap <- umap
# DEA --------------------------------------------------------------------------
myscobj <- subset(scobj, subset = disease__ontology_label %in% c("chronic obstructive pulmonary disease", "normal"))
myscobj <- subset(myscobj, subset = cell_type__ontology_label != "apoptosis fated cell")
myscobj <- subset(myscobj, subset = cell_type__ontology_label != "mitotic cell cycle")
myscobj@assays$RNA$data = myscobj@assays$RNA$counts
myscobj$health <- factor(myscobj$health, levels = c("healthy","COPD"))

# table(myscobj$health)
# COPD healthy 
# 68792  259912 
# cell - cell ==================================================================
all.markers <- FindAllMarkers(myscobj, group.by = "cell_type__ontology_label")

save(all.markers, file = "findallmarkers.Rdata")
write.csv(all.markers, file = "findallmarkers.csv")
# health - COPD ================================================================
suppressWarnings(require(crayon, warn.conflicts = F))

Pa    = crayon::cyan
Er    = crayon::red$bold
Sa    = crayon::blue
No    = crayon::magenta$bold
Wa    = crayon::yellow

mymeta <- myscobj@meta.data %>% 
  dplyr::select(biosample_id, donor_id, health) %>% 
  distinct(biosample_id, donor_id, health, .keep_all = TRUE) %>% 
  mutate(biosample_id = str_replace_all(biosample_id, "_", "-"))

for (i in unique(myscobj$cell_type__ontology_label)) {
  tmp <- myscobj %>% subset(subset = cell_type__ontology_label == i)
  
  counts   <- AggregateExpression(tmp, group.by = "biosample_id")
  counts   <- counts$RNA %>% as.data.frame
  colnames(counts) <- str_remove(colnames(counts), "^g")
  mymetacon<- mymeta %>% filter(biosample_id %in% colnames(counts))
  counts   <- counts[,mymetacon$biosample_id]
  condition<- factor(mymetacon$health, levels = c("healthy","COPD"))
  
  message(Sa(paste0("--> Start ", i," DEA <--")))
  if(all(mymetacon$biosample_id == colnames(counts))){
    DEG <- data.frame(
      symbol = rownames(counts),
      health = rowMeans(counts[,condition == "healthy"]),
      copd   = rowMeans(counts[,condition == "COPD"]),
      t      = apply(counts, 1, function(x){t.test(x ~ condition, var.equal = TRUE)$p.value}),
      wilcox = apply(counts, 1, function(x){wilcox.test(x ~ condition)$p.value}),
      stringsAsFactors = FALSE
    ) %>% 
      mutate(FC = copd / health, 
             log2FC = log2(copd / health),
             fdrt = p.adjust(t, method = "fdr"),
             fdrwilcox = p.adjust(wilcox, method = "fdr"))
    
    save(counts, file = paste0(i,"_counts.Rdata"))
    save(DEG,    file = paste0(i,"_DEGs.Rdata"))
    write.csv(DEG, paste0(i,"_DEGs.csv"))
  }
  message(Sa(paste0("--> End ", i," DEA <--")))
}

pList <- list()
n <- 1
for (i in c("B cell","Endothelial","Epithelial","Mast","Mesenchyme","Myeloid",
            "Neutrophil","Plasma","T cell")) {
  load(paste0(i,"_DESeq2-DEGs.Rdata"))
  p <- geom_volcano(DEG_drinker_vs_nondrinker, title = i)
  pList[[n]] <- p
  n <- n + 1
}

p_out <- patchwork::wrap_plots(pList) +
  plot_layout(ncol = 3, nrow = 3)
ggsave(p_out, filename = paste0(i,"volcanoplot.pdf"), width=15, height=15, units="in")
# Visualization ----------------------------------------------------------------
table(myscobj$cell_type__ontology_label)
# B cell 
# 6769 
# endothelial cell 
# 11650 
# epithelial cell 
# 6322 
# fibroblast 
# 15873 
# lymph node lymphatic vessel endothelial cell 
# 2701 
# lymphocyte 
# 46887 
# mast cell 
# 2898 
# myeloid cell 
# 180703 
# smooth muscle cell 
# 4995 
# type I pneumocyte 
# 18855 
# type II pneumocyte 
# 31051

markgenes <- c("CD79A","VWF","EPCAM","COL1A1","DCN","LYVE1","CD3E","NKG7","TPSAB1",
               "CD14","FABP4","RGS5","ACTA2","SFTPC","AGER")
focusgene <- c("SUCLG1","SUCLA2","SDHA","SDHB","SDHC","SDHD","OGDH",
               "ALDH5A1","ALAS1","ALAS2","SUCNR1") # "SCULG2", "MUT",
# UMAP ----
p1 <- DimPlot(myscobj, reduction = "umap", group.by = "cell_type__ontology_label", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(myscobj, reduction = "umap", group.by = "health", cols = c("#599CB4", "#C25759"), 
              label = TRUE, pt.size = 0.1, raster = FALSE)
p <- p1 + p2 + plot_layout(ncol = 1, nrow = 2)
ggsave(p, filename = "UMAP.pdf", width=8, height=10, units="in")
# Alluvium ----
p <- myscobj@meta.data %>% 
  mutate(value = 1) %>% 
  group_by(health, cell_type__ontology_label) %>%
  summarise(n = n()) %>% 
  group_by(health) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = health, y = n, 
             fill = cell_type__ontology_label, 
             stratum = cell_type__ontology_label, 
             alluvium = cell_type__ontology_label))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Alluvium.pdf",
       width=10, height=10, units="in")

# FeaturePlot ----
p_vln <- VlnPlot(myscobj, features = focusgene, cols = colors_list, group.by = "cell_type__ontology_label")
ggsave(p_vln, filename = "VlnPlot.pdf", width=16, height=16, units="in")

p_vln_epi <- VlnPlot(subset(myscobj, subset = cell_type__ontology_label == "epithelial cell"), 
                     features = focusgene, cols = colors_list, group.by = "health")
ggsave(p_vln_epi, filename = "VlnPlot(epi).pdf", width=12, height=16, units="in")
p_vln_fibro<- VlnPlot(subset(myscobj, subset = cell_type__ontology_label == "fibroblast"), 
                     features = focusgene, cols = colors_list, group.by = "health")
ggsave(p_vln_fibro, filename = "VlnPlot(fibro).pdf", width=12, height=16, units="in")

p_out <- lapply(focusgene, 
                function(x){FeaturePlot(myscobj, features = x) + 
                    scale_colour_gradientn(
                      colours = colorRampPalette(
                        c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                          '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                          '#6F003D','#56033F'))(1000))}) %>% 
  patchwork::wrap_plots(ncol = 3, nrow = 4)
ggsave(p_out, filename = "FeaturePlot.pdf", width=14, height=16, units="in")

# DotPlot ----
p <- DotPlot(myscobj, features = markgenes,
  group.by = "cell_type__ontology_label") + 
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "MARKER.pdf", width=6, height=8, units="in")
# volcano ----
geom_volcano <- function(dat,pos.num = 10,neg.num = -10, title){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    # tibble::rownames_to_column("SYMBOL") %>% 
    filter(complete.cases(log2FC)) %>% 
    mutate(color = case_when(log2FC >  1 & fdrt < 0.05 ~ 2,
                             log2FC < -1 & fdrt < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FC))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$symbol[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c("",0,0,1,1,1,0,1,1,NA,NA)
    dat_up$log2FC <- dat_up$log2FC %>% as.numeric
    dat_up$fdrt <- dat_up$fdrt %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FC)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$symbol[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,1,1,1,0,1,1,NA,NA)
    dat_down$log2FC <- dat_down$log2FC %>% as.numeric
    dat_down$fdrt <- dat_down$fdrt %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FC, y = -log10(fdrt), col = log2FC, label = symbol)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
    )+
    geom_vline(xintercept = c(-1,1), color = "gray80", linetype = 2) +
    geom_hline(yintercept = c(1.30103), color = "gray80", linetype = 2) +
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab("Log2(Fold Change)") +
    labs(color = "") +
    scale_color_gradient2(high = "red3", mid = "white", low = "blue3",
                          midpoint = 0, na.value = "grey80"
                          #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
    )+
    scale_size_continuous(range = c(0.1, 4)) +
    geom_text_repel(
      data = dat_up,
      color = "red3",
      size = 5,
      nudge_x = pos.num - as.numeric(dat_up$log2FC),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FC),
      segment.size = 0.3,
      segment.color = "grey",
      direction="y",
      hjust= 1,
      max.overlaps = Inf) +
    labs(title = title) + 
    # labs(size = expression("Abundance (log2)"),
    #      color = expression("Direction signed"),
    #      title = trait.names[i]) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
          #legend.spacing.y = unit(1,"cm")
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, size=20))
}

pList <- list()
n <- 1
for (i in c("B cell",                                      
            "endothelial cell",                            
            "epithelial cell",                             
            "fibroblast",                                  
            "lymph node lymphatic vessel endothelial cell",
            "lymphocyte",                                  
            "mast cell",                                   
            "myeloid cell",                                
            "smooth muscle cell",                          
            "type I pneumocyte",                           
            "type II pneumocyte")) {
  load(paste0(i,"_DEGs.Rdata"))
  p <- geom_volcano(DEG, title = i)
  pList[[n]] <- p
  n <- n + 1
}

p_out <- patchwork::wrap_plots(pList) +
  plot_layout(ncol = 3, nrow = 6)
ggsave(p_out, filename = "volcanoplot.pdf", width=14, height=24, units="in")
# correlation ----
hallmarks <- list(
   qusage::read.gmt("FOROUTAN_INTEGRATED_TGFB_EMT_DN.v2025.1.Hs.gmt"),
   qusage::read.gmt("FOROUTAN_INTEGRATED_TGFB_EMT_UP.v2025.1.Hs.gmt"),
   qusage::read.gmt("FOROUTAN_TGFB_EMT_DN.v2025.1.Hs.gmt"),
   qusage::read.gmt("FOROUTAN_TGFB_EMT_UP.v2025.1.Hs.gmt")
)
names(hallmarks) <- c("INT_TGFB_EMT_DN","INT_TGFB_EMT_UP","TGFB_EMT_DN","TGFB_EMT_UP")

myscobj <- AddModuleScore(
  object = myscobj, 
  features = hallmarks, 
  name = c("INT_TGFB_EMT_DN","INT_TGFB_EMT_UP","TGFB_EMT_DN","TGFB_EMT_UP"),     
  seed = 1011,                
  ctrl = 100                
)
emt_coranalysis <- lapply(focusgene, function(gene){
  lapply(c("INT_TGFB_EMT_DN1","INT_TGFB_EMT_UP2","TGFB_EMT_DN3","TGFB_EMT_UP4"), 
         function(pathway){
           res <- cor.test(myscobj[[pathway]] %>% unlist, 
                           myscobj@assays$RNA$data[gene,] %>% unlist)
           data.frame(
             gene    = gene,
             pathway = pathway,
             est     = res$estimate,
             pval    = res$p.value
           )
         }) %>% bind_rows()
}) %>% bind_rows()
























