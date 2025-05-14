require(corrplot)

dat <- readxl::read_excel("41.correlation heatmap.xlsx") %>% 
  tibble::column_to_rownames("Sample")

dat.clust <- dat %>% t %>% 
  dist() %>% 
  hclust(method = "ward.D2")
dat <- dat[,dat.clust$order]

mycor <- dat %>% cor
corrplot(mycor,
         tl.col = "#222222", tl.pos = "td",
         type = "upper", 
         pch.cex = 0.5,
         p.mat = mycor,
         diag = F, 
         method = "color")