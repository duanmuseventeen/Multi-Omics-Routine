# Ref:
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# https://github.com/ixxmu/mp_duty/issues/3901
# https://www.jianshu.com/p/e85d23a25a43
# https://www.jianshu.com/p/e30c307473a8

# 软件整理
# 1. RNA velocity
# 2. Monocle3 (Monocle2)
# 3. scVelo
# 4. slingshot  
# 5. velocyto
# 6. Dynamo
# 7. tradeSeq

scobj.sce <- as.SingleCellExperiment(scobj, assay = "RNA")

# Method 1----------------------------------------------------------------------
scobj.sim <- slingshot(scobj.sce, clusterLabels = 'cell_type', reducedDim = 'UMAP',
                       start.clus= "Tn", end.clus = NULL)

# scobj.sce@colData %>% View
# summary(scobj.sce$slingPseudotime_1)

if(FALSE){
  # 保存时，把下面的注释运行
  # pdf(file = "scobj lineage.pdf",width = 7, height = 5)
  # colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  # 
  # plotcol <- colors[cut(scobj.sim$slingPseudotime_1, breaks=100)]
  # plot(reducedDims(scobj.sim)$UMAP, col = plotcol, pch=16, asp = 1)
  # lines(SlingshotDataSet(scobj.sim), lwd=2, col=brewer.pal(9,"Set1"))
  # legend("right",
  #        legend = paste0("lineage",1:6),
  #        col = unique(brewer.pal(6,"Set1")),
  #        inset=0.8,
  #        pch = 16)
  # 
  # plotcol <- colors[cut(scobj.sim$slingPseudotime_2, breaks=100)]
  # plot(reducedDims(scobj.sim)$UMAP, col = plotcol, pch=16, asp = 1)
  # lines(SlingshotDataSet(scobj.sim), lwd=2, col=brewer.pal(9,"Set1"))
  # legend("right",
  #        legend = paste0("lineage",1:6),
  #        col = unique(brewer.pal(6,"Set1")),
  #        inset=0.8,
  #        pch = 16)
  # dev.off()
}

# Graphics with ggplot2 system
dat <- reducedDims(scobj.sce)$UMAP %>% as.data.frame %>%
  tibble::rownames_to_column("barcode") %>% as.data.frame %>%
  left_join(fig@meta.data %>%
              dplyr::select("cell_type") %>%
              tibble::rownames_to_column("barcode"),
            by = "barcode")
dat.fit1 <- scobj.sim@colData@listData$slingshot@metadata$curves$Lineage1
dat.fit2 <- scobj.sim@colData@listData$slingshot@metadata$curves$Lineage2

dat$slingPseudotime_1 <- scobj.sim$slingPseudotime_1
dat$slingPseudotime_2 <- scobj.sim$slingPseudotime_2

# p <- ggplot() +
#   geom_point(data = dat %>% filter(complete.cases(slingPseudotime_1)),
#              aes(x = umap_1, y = umap_2, col = slingPseudotime_1, alpha = 0.1)) +
#   scale_color_viridis(name="Slingshot:\nPseudotime1", option = "D") +
#   ggnewscale::new_scale_color() +
#   geom_point(data = dat %>% filter(complete.cases(slingPseudotime_2)),
#              aes(x = umap_1, y = umap_2, col = slingPseudotime_2, alpha = 0.1)) +
#   geom_path(data = dat.fit1$s[dat.fit1$ord,],
#             aes(x = umap_1, y = umap_2), col = "red3") +
#   geom_path(data = dat.fit2$s[dat.fit2$ord,],
#             aes(x = umap_1, y = umap_2), col = "blue3") +
#   scale_color_viridis(name="Slingshot:\nPseudotime2", option = "A") +
#   guides(colour = guide_legend(title = "Slingshot:\nPseudotime2")) +
#   ggtitle("Slingshot scobj") +
#   theme_bw() +
#   theme(text = element_text(size = 16),
#         plot.margin = ggplot2::margin(15,15,15,15),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5))
# ggsave(filename = "scobj cell(Slingshot).pdf", plot = p, device = "pdf",
#        width = 7, height = 5, dpi = 300, units = "in")

p1 <- ggplot() +
  geom_point(data = dat %>% filter(complete.cases(slingPseudotime_1)),
             aes(x = umap_1, y = umap_2, col = slingPseudotime_1)) +
  scale_color_viridis(name="Slingshot:\nPseudotime1", option = "A") +
  geom_path(data = dat.fit1$s[dat.fit1$ord,],
            aes(x = umap_1, y = umap_2), col = "red3") +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(-8,8)) +
  ggtitle("Slingshot scobj") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = ggplot2::margin(15,15,15,15),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
p2 <- ggplot() +
  geom_point(data = dat %>% filter(complete.cases(slingPseudotime_2)),
             aes(x = umap_1, y = umap_2, col = slingPseudotime_2)) +
  scale_color_viridis(name="Slingshot:\nPseudotime2", option = "D") +
  geom_path(data = dat.fit2$s[dat.fit1$ord,],
            aes(x = umap_1, y = umap_2), col = "green4") +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(-8,8)) +
  ggtitle("Slingshot scobj") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = ggplot2::margin(15,15,15,15),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
p_out <- p1 + p2
ggsave(filename = "scobj cell(Slingshot).pdf", plot = p_out, device = "pdf",
       width = 12, height = 5, dpi = 300, units = "in")

# Method 2: step by step--------------------------------------------------------
rd <- dat %>% tibble::column_to_rownames("barcode") %>% dplyr::select(umap_1, umap_2)
pto <- getLineages(rd, dat$cell_type, start.clus = 'Tn_CCR7')
plot(rd, col = brewer.pal(9,"Set1")[
  as.numeric(factor(dat$cell_type, 
                    levels = c("Tn", "CD8 T unknown","CD8 Teff",
                               "CD8 TemRA", "CD8 Tex 1","CD8 Tex 2",
                               "CD8 Trm")))
], asp = 1, pch = 16)
lines(SlingshotDataSet(pto), lwd = 3, col = 'black')

crv1 <- getCurves(pto)
crv1
plot(rd, col = brewer.pal(9,"Set1")[
  as.numeric(factor(dat$cell_type, levels = c("Tn_CCR7", "CD8 T unknown_AOAH","CD8 Teff_XCL1/2",
                                              "CD8 TemRA_FGFBP2", "CD8 Tex_HAVCR2","CD8 Tex_SPRY1",
                                              "CD8 Trm_NR4A2/3")))
],
asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')