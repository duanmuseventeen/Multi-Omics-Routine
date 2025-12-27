# data is a numeric df with rowname, in which row means sample and column means feature
# group should be factor

pca <- function(dat, group){
  
  suppressMessages(require(ggplot2))
  suppressMessages(require(FactoMineR))
  
  pre.pca <- PCA(dat, graph = FALSE)
  factoextra::fviz_pca_ind(pre.pca,
                           mean.point = F,
                           pointsize = 4,
                           geom= "point",
                           col.ind = group,
                           addEllipses = F) +
    stat_ellipse(level = 0.95, type = "norm", geom = "polygon", col = "gray20", alpha = 0) +
    # geom_point(size = 4) +
    # ggtitle("") +
    # scale_shape_manual(values = ) +
    # scale_fill_manual(values = ) +
    # scale_color_manual(values = ) +
    guides(size = "none") +
    labs(col = "Group", shape = "Group") +
    ggtitle("") +
    theme_bw() +
    theme(text = element_text(size = 20),
          plot.margin =  ggplot2::margin(30,30,30,30),
          panel.grid = element_blank(),
          legend.background = element_rect(linetype = 1, color = "black"),
          legend.position = "right",
          axis.text = element_text(size= 16, color = "black"))
}
