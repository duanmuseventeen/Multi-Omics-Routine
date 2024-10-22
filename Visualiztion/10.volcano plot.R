# plot volcano
# Not run

# Proteomic predictors of individualized nutrient-specific insulin secretion in health and disease
# PMID: 38959864

require(ggplot2)
require(ggrepel)

# dat is a df including log2FC, padj, color (equal to log2FC but NA when padj > 0.05)
# and abundance
geom_volcano <- function(dat){
  ggplot(dat, aes(x = log2FC, y = -log10(padj), col = color, label = label)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
      )+
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab(expression("R")) +
    scale_color_gradient2(high = "red3", mid = "white", low = "blue3",
                          midpoint = 0, na.value = "grey80"
                          #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
    )+
    scale_size_continuous(range = c(0.1, 4))+
    geom_text_repel(
      data = subset(dat, log2FC > 0), 
      color = "red",
      size = 5,
      nudge_x = 0.02 - subset(dat, log2FC > 0)$log2FC,
      segment.size=0.3, 
      segment.color="grey",
      direction="y", 
      hjust= 0,
      max.overlaps = Inf)+
    geom_text_repel(
      data=subset(dat, log2FC < 0), 
      color="blue",
      size=5,
      nudge_x = -0.02 - subset(dat, log2FC < 0)$R,
      segment.size = 0.3, 
      segment.color = "grey", 
      direction="y", 
      hjust= 1,
      max.overlaps = Inf)+
    labs(size = expression("Abundance (log2)"),
         color = expression("Direction signed"),
         title = trait.names[i]) +  
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
          #legend.spacing.y = unit(1,"cm")
    ) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "italic", colour="grey50", size=20))
}

ggsave(filename = "volcanoplot.pdf", width=10, height=6.5, units="in")