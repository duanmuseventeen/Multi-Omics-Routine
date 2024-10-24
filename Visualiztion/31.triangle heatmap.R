# triangle heatmap

# dat with rowname is a df, in which row means sample and column means feature
dat <- matrix(runif(200) , nrow = 10) %>% 
  as.data.frame

Features <- colnames(dat)

datcor <- cor(dat) %>% as.data.frame
datcorsig <- cor.mtest(dat, conf.level = 0.95, method = "pearson")$p
datcorsig.fdr <- datcorsig %>% 
  p.adjust(method = "fdr") %>%
  matrix(nrow = length(Features)) %>% 
  as.data.frame 

datcor[row(datcor) >= col(datcor)] <- NA
datcorsig.fdr[row(datcorsig.fdr) >= col(datcorsig.fdr)] <- NA

est <- datcor %>% 
  as.data.frame %>% 
  tibble::rownames_to_column("A") %>% 
  tidyr::pivot_longer(cols = Features, 
                      names_to = "B", values_to = "cor") %>% 
  filter(!is.na(cor)) %>%
  mutate(AB = paste0(A, "-", B),
         A = factor(A, levels = Features, labels = Features),
         B = factor(B, levels = Features, labels = Features)) 

rownames(datcorsig.fdr) <- colnames(datcorsig.fdr) <- colnames(datcorsig)
sig <- datcorsig.fdr %>% 
  as.data.frame %>% 
  tibble::rownames_to_column("A") %>% 
  tidyr::pivot_longer(cols = Features, 
                      names_to = "B", values_to = "fdr") %>% 
  filter(!is.na(fdr)) %>%
  # filter(fdr < 0.05) %>% 
  mutate(AB = paste0(A, "-", B),
         A = factor(A, levels = Features, labels = Features),
         B = factor(B, levels = Features, labels = Features)) %>% 
  dplyr::select(-c(A,B))

p <- 
  est %>% 
  left_join(sig, by = "AB") %>% 
  mutate(lab = case_when(fdr > 0.05 ~ "",
                         fdr > 0.01 ~ "*",
                         fdr > 0.001 ~ "**",
                         TRUE ~ "***")) %>% 
  ggplot(aes(A, B, label = lab)) +
  geom_tile(aes(fill = cor)) +
  geom_text(angle = 45) +
  coord_equal(clip = "off") +
  scale_fill_gradient2(
    low = "blue3",
    mid = "white",
    high = "red3",
    # na.value = "white",
    midpoint = 0
  ) + 
  labs(x = "", y = "") +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 1), 
        plot.margin = ggplot2::margin(80,80,80,80),
        legend.position = "top")

p.rot45 <- ggplotify::as.ggplot(p, angle = -45)
p.rot45
