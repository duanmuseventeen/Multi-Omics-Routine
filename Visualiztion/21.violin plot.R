# Violin plot

require(ggplot2)
require(dplyr)

dat <- readxl::read_excel("8-38347143-Source Data Extended Fig8.xlsx", sheet = "ED_8b") %>% 
  mutate(x = factor(Dimension, 
                    levels = c("Single","Mid","Multi"), 
                    labels = c("1", "2-3", "4-6"))) %>% 
  ggplot(aes(x = x, y = `C-index`, fill = x)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1) +
    annotate(geom = "text", 
             x = c(1.5, 2.5), 
             y = c(0.84,0.88), 
             label = c("FDR = 0.01", "FDR = 0.11"), size = 6) +
    annotate(geom = "segment", 
             x = c(1,2), xend = c(2,3), 
             y = c(0.82,0.86), yend = c(0.82,0.86)) +
    annotate(geom = "text",
             x = 2, y = 0.9,
             parse = TRUE, 
             label = "italic(P) == 0.0027", size = 6) +
    scale_fill_manual(values = c("#dfe3ec","#a3b2c9","#6683a7")) +
    labs(x = "Amount of data dimension", y = "C index in test set") +
    theme_classic() +
    theme(text = element_text(size = 16), legend.title = element_blank(),
          plot.margin = ggplot2::margin(30,30,30,30),
          panel.grid = element_blank(), legend.background = element_blank(),
          axis.title.x = element_text(vjust = -4),
          axis.title.y = element_text(vjust = 4))
