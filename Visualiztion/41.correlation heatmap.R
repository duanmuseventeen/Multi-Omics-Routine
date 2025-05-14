require(dplyr)
require(ggplot2)

gut <- readxl::read_excel("41.correlation heatmap.xlsx") %>% 
  tibble::column_to_rownames("Sample")

ggcorheat <- function(dat, cor.method = "spearman", 
                      colours = c("blue3","white","red3")){
  cor_est_dat <- cor_sig_dat <- matrix(ncol = ncol(dat), nrow = ncol(dat))
  colnames(cor_sig_dat) <- rownames(cor_sig_dat) <-
    colnames(cor_est_dat) <- rownames(cor_est_dat) <- 
    colnames(dat)
  
  for (i in colnames(dat)) {
    for (j in colnames(dat)) {
      tmp <- cor.test(dat %>% dplyr::select(i) %>% unlist, 
                      dat %>% dplyr::select(j) %>% unlist,
                      method = cor.method)
      cor_sig_dat[j,i] <- tmp$p.value
      cor_est_dat[j,i] <- tmp$estimate
    }
  }
  cor_sig_dat <- cor_sig_dat %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("ID") %>% 
    tidyr::pivot_longer(cols = colnames(dat)) %>% 
    mutate(padj = p.adjust(value, method = "fdr")) %>% 
    dplyr::select(-value) %>% 
    tidyr::pivot_wider(names_from = "name", values_from = "padj") %>% 
    tibble::column_to_rownames("ID")
  
  datagg.p <- cor_sig_dat %>% as.data.frame %>% 
    tibble::rownames_to_column("Gene") %>%
    reshape2::melt("Gene") %>%
    rename(p = "value") %>%
    mutate(ID = paste0(Gene,"-",variable)) %>%
    dplyr::select(-c(Gene, variable)) %>%
    mutate(sig = ifelse(p > 0.05,"",ifelse(p > 0.01, "*",ifelse(p > 0.001, "**", "***"))))
  datagg <- cor_est_dat %>% as.data.frame %>% 
    tibble::rownames_to_column("Gene") %>%
    reshape2::melt("Gene") %>%
    mutate(ID = paste0(Gene,"-",variable)) %>%
    left_join(datagg.p, by = "ID")
  
  data.clust <- datagg %>% 
    dplyr::select(all_of(c("Gene","variable","value"))) %>% 
    tidyr::pivot_wider(names_from = "variable", values_from = "value") %>% 
    tibble::column_to_rownames("Gene") %>% 
    dist() %>% 
    hclust(method = "ward.D2")
  
  p <- datagg %>% 
    mutate(Gene = factor(Gene, levels = colnames(dat)[data.clust$order]),
           variable = factor(variable, levels = colnames(dat)[data.clust$order])) %>% 
    ggplot(aes(x = Gene, y = variable, fill = value, label = sig)) +
    geom_tile() +
    geom_text(color = "white") +
    labs(x = "", y = "", fill = "spearman") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient2(low = colours[1], mid = colours[2], high = colours[3], 
                         midpoint = 0, breaks = c(-1,-.5,0,.5,1), limits = c(-1,1)) +
    coord_fixed(ratio = 1) + 
    theme_classic() +
    theme(text = element_text(size = 12),
          axis.text.x.bottom = element_text(angle = 270, hjust = 0),
          line = element_blank())
  return(p)
}

p1 <- ggcorheat(dat = gut)
ggsave(plot = p1, filename = paste0("gut",".pdf"), 
       device = "pdf",width = 8, height = 8, dpi = 300)

# Interestingly, even we conduct FDR to avoid the false positive results, particularly,
# in correlation analysis, FDR adjustion increase the number of significant correlation (
# fdr-q < 0.05).



