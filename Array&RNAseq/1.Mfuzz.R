require(Mfuzz)
require(RColorBrewer)

set.seed(1011)
dat <- # counts, fpkm, tpm or samilar things
n.group <- 
  
eset <- matrix(nrow = nrow(fpkm), ncol = n.group)
eset[,1] <- dat %>% dplyr::select(sample.names.in.group1) %>% rowMeans
eset[,2] <- dat %>% dplyr::select(sample.names.in.group2) %>% rowMeans
eset[,3] <- dat %>% dplyr::select(sample.names.in.group3) %>% rowMeans
eset[,4] <- dat %>% dplyr::select(sample.names.in.group4) %>% rowMeans
...

colnames(eset) <- group.names
rownames(eset) <- rownames(dat)
eset <- new("ExpressionSet", exprs = eset)
eset <- standardise(eset)

c <- 9 
m <- mestimate(eset) 
cl <- mfuzz(eset, c = c, m = m) # 聚类
cl$size
cl$membership

Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
mfuzz.plot(
  eset,
  cl,
  mfrow=c(3,3),
  new.window= FALSE,
  time.labels=colnames(eset)
  )

mfuzz.plot2(
  eset,
  cl,
  mfrow=c(2,3),
  time.labels=colnames(eset),
  ylim.set=c(0,0), 
  xlab = "",
  ylab = "Expression changes",
  x11 = FALSE #If TRUE, a new window will be open for plotting.
)

dMfuzz <- data.frame(
  `Gene Symbol` = names(cl$cluster),
  cluster = cl$cluster,
  row.names = names(cl$cluster),
  stringsAsFactors = F
)