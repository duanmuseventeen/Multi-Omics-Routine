# volcano plot

# for one dataset
geom_volcano <- function(dat){
  ggplot(dat, aes(x = logFC, y = -log10(adj.P.Val), col = sig)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, col = "gray80") +
    # geom_vline(xintercept = c(-1, 1), linetype = 2, col = "gray80") +
    geom_point() +
    scale_color_manual(values = c("gray80","#0C56A0","#b9292b")) + 
    ylab(expression("-log10("~italic(P~value)~")")) +
    guides(color = "none") +
    theme_classic() +
    theme(text = element_text(size = 20), 
          panel.grid = element_blank(), 
          plot.margin = ggplot2::margin(40,40,40,40), ,
          legend.background = element_rect(linetype = 1, color = "black"),
          legend.position = c(0.8,0.8),
          axis.text = element_text(size= 16, color = "black"))
}

# for multi datasets
ggplot(dat, aes(x = dataset_ID, y = log2FC, colour = sig)) +
  geom_jitter() +
  scale_color_manual(values = c("gray80","#0C56A0","#b9292b")) +
  labs(x = "", y = "Log2(FC)") + 
  theme_bw() +
  theme(text = element_text(size = 20), 
        panel.grid = element_blank(), 
        plot.margin = ggplot2::margin(40,40,40,40), ,
        legend.background = element_rect(linetype = 1, color = "black"),
        legend.position = c(0.8,0.8),
        axis.text = element_text(size= 16, color = "black"))