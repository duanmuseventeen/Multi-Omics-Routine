require(corrplot)

dat <- readxl::read_excel("41.correlation heatmap.xlsx") %>% 
  tibble::column_to_rownames("Sample")

cor.method = "spearman"

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
# cor_sig_dat <- cor_sig_dat %>% 
#   as.data.frame %>% 
#   tibble::rownames_to_column("ID") %>% 
#   tidyr::pivot_longer(cols = colnames(dat)) %>% 
#   mutate(padj = p.adjust(value, method = "fdr")) %>% 
#   dplyr::select(-value) %>% 
#   tidyr::pivot_wider(names_from = "name", values_from = "padj") %>% 
#   tibble::column_to_rownames("ID")

corrplot(cor_est_dat %>% as.matrix,
         tl.col = "#222222", 
         tl.pos = "td",
         type = "upper",
         pch.cex = 0.5,
         p.mat = cor_sig_dat %>% as.matrix,
         sig.level = 0.05,
         diag = F,
         method = "color")