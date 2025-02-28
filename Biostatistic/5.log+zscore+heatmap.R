# http://www.ehbio.com/Bioinfo_R_course/Rplots.html#ggplot2_heatmapnonlinear
# https://zhuanlan.zhihu.com/p/576687876
# https://cloud.tencent.com/developer/article/1486102
# https://blog.csdn.net/leianuo123/article/details/107669123

# what fun scale do?------------------------------------------------------------
a <- runif(100)

a.scale <- scale(a)
a.scale

zscore <- (a - mean(a)) / sd(a)
zscore

a.scale - zscore

# deal with outliers in heatmap-------------------------------------------------
# the data has been log2 transformed

load("5.log+zscore+heatmap.Rdata")

dim(dat)
# [1] 127 168
str(dat)

dat <- dat %>% 
  dplyr::select(-ions) %>% 
  tibble::column_to_rownames("name")
dat.num <- apply(dat, 2, as.numeric)
dat.num <- as.data.frame(dat.num)
rownames(dat.num) <- rownames(dat)

dim(dat.num)
# [1] 127 168
str(dat.num)

dat.scale <- apply(dat.num, 1, scale)
dat.scale <- dat.scale %>% t %>% as.data.frame
colnames(dat.scale) <- colnames(dat)

dim(dat.scale)
# [1] 127 168
str(dat.scale)

pheatmap::pheatmap(dat.scale, show_colnames = FALSE, scale = "none")

summary(unlist(dat.scale))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -5.56638 -0.64017  0.06696  0.00000  0.64531  7.06820 

# delete 1% and 99%-------------------------------------------------------------
quantile(unlist(dat.scale), c(.01,.99))

dat.adj <- dat.scale
dat.adj[dat.adj > 2.270816] <- 2.270816
dat.adj[dat.adj < -2.531722]<- -2.531722

pheatmap::pheatmap(dat.adj, clustering_method = "ward.D2",
                   show_colnames = FALSE, 
                   scale = "none")
# delete 2% and 98%-------------------------------------------------------------
quantile(unlist(dat.scale), c(.02,.98))

dat.adj <- dat.scale
dat.adj[dat.adj > 1.973409] <- 1.973409
dat.adj[dat.adj < -2.146326]<- -2.146326

pheatmap::pheatmap(dat.adj, clustering_method = "ward.D2",
                   show_colnames = FALSE, 
                   scale = "none")

# delete 5% and 95%-------------------------------------------------------------
quantile(unlist(dat.scale), c(.05,.95))

dat.adj <- dat.scale
dat.adj[dat.adj > 1.527590] <- 1.527590
dat.adj[dat.adj < -1.701551]<- -1.701551

pheatmap::pheatmap(dat.adj, clustering_method = "ward.D2",
                   show_colnames = FALSE, 
                   scale = "none")
