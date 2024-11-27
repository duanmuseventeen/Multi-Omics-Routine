require(mixOmics)

# dat.dea is numeric df, in which row means sample and column means featrue
# group should be factor

dat.pls <- plsda(dat.dea, group$group, ncomp = 3, scale = TRUE)

plotIndiv(dat.pls, ind.names = TRUE, ellipse = TRUE, legend = TRUE)

ggplot(data = dat.pls$variates$X,
       aes(x = comp1,y = comp2, colour = group, fill = group)) +
  geom_hline(yintercept = 0, color = "#999999", linetype = 2) +
  geom_vline(xintercept = 0, color = "#999999", linetype = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1) +
  geom_point(size = 1, alpha = 1)
  # scale_color_manual() +
  guides(fill = "none") +
  labs(x = paste0("PC1: ", sprintf("%.2f", dat.pls$prop_expl_var$X[1] * 100),"%"),
       y = paste0("PC2: ", sprintf("%.2f", dat.pls$prop_expl_var$X[2] * 100),"%")) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.title = element_blank(),
        plot.margin = ggplot2::margin(30,30,30,30),
        panel.grid = element_blank(), legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = c(0.2,0.8))

ggsave(filename = "PLS.pdf", path = "D:/", device = "pdf", width = 5, height = 5, dpi = 300)
