# optimization of scRNAseq graphics
# use data from GSE115469

# https://mp.weixin.qq.com/s/CljLXy9au2sjm4sR_gp9wQ

# R-----------------------------------------------------------------------------
## load pkgs---------------------------------------------------------------------
library(tidydr) # Y叔
library(ggplot2)
library(scCustomize)
## 1 axis------------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/669051636
DimPlot(scobj.h.sc, reduction = "UMAP", group.by = "cell_type", label = TRUE, pt.size = 1.2) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

## 2 highlight cluster-----------------------------------------------------------
# method 1:
# https://blog.csdn.net/qq_42090739/article/details/123470793
# please do that in AI or PowerPoint

# method 2:
# https://github.com/sajuukLyu/ggunchull
# devtools::install_github("sajuukLyu/ggunchull", type = "source")

require(ggunchull)

dat.plot <- scobj.h.sc@reductions$UMAP@cell.embeddings %>% 
  as.data.frame %>% 
  mutate(cell_type = scobj.h.sc@meta.data$cell_type)

delta_r <- 0.03
th_r <- 0.03
delta <- diff(range(dat.plot$umap_1)) * delta_r
th <- diff(range(dat.plot$umap_1)) * th_r

ggplot(dat.plot, aes(x = umap_1, y = umap_2, fill = cell_type, color = cell_type)) +
  stat_unchull(
    # data = subset(dat.plot, celltype == "Hepatocyte"),
    alpha = 0.2, size = 0.5,lty = 2, delta = delta, th = th) +
  geom_point(size = 1.5,show.legend = F) +
  guides(fill = guide_legend(override.aes = list(colour = "white", linetype = 1, size = 1))) +
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = 2,hjust = 0.03))
# method 3----------------------------------------------------------------------
DimPlot(scobj.h.sc, reduction = "UMAP", group.by = "cell_type", label = TRUE, 
        pt.size = 1.2, 
        cols.highlight = "pink",
        cells.highlight = rownames(scobj.h.sc@meta.data)[scobj.h.sc@meta.data$cell_type == "Hepatocyte"]) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
## 3 label-----------------------------------------------------------------------
require(ggrepel)

dat.plot <- scobj.h.sc@reductions$UMAP@cell.embeddings %>% 
  as.data.frame %>% 
  mutate(cell_type = scobj.h.sc@meta.data$cell_type)

cell_type.pos <- dat.plot %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
    )

DimPlot(scobj.h.sc, reduction = "UMAP", group.by = "cell_type", label = FALSE, pt.size = 1.2) + 
  geom_label_repel(aes(x = umap_1,y = umap_2,label = cell_type, color = cell_type), 
                   fontface = "bold",
                   data = cell_type.pos, # locate labels
                   box.padding = 0.5) +
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() 

## 4 Contours-------------------------------------------------------------------
# Refer to 39346916, in my opinion, contours match t-sne better
DimPlot(scobj.h.sc, reduction = "TSNE", group.by = "cell_type", pt.size = 1.5)+
  stat_density_2d(aes(x = tSNE_1, y = tSNE_2), linemitre = 10, col = "black")
## 5 Density--------------------------------------------------------------------
# BiocManager::install("Nebulosa")
require(Nebulosa)

plot_density(scobj.h.sc, "ACTB")
plot_density(scobj.h.sc, c("ACTB","GAPDH"), joint = TRUE)
## *example----------------------------------------------------------------------
# https://cloud.tencent.com/developer/article/2317517
# python------------------------------------------------------------------------
# https://blog.csdn.net/wcy1995427/article/details/140009111
# https://www.cnblogs.com/starlitnightly/p/18262565
