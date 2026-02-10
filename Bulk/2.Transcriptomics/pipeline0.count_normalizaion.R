# Ref: 
# https://mp.weixin.qq.com/s/imJ5za173pKFJVad5j_5iw
# https://www.jianshu.com/p/73cb2ebcd1b0

countToTpm <- function(counts, effLen){
  rate <- log(counts) - log(effLen) # row respresents genes, align with effLen
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


library(GenomicFeatures)
## 1.读取GTF文件
txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf")
## 2.提取外显子信息
exonic <- exonsBy(txdb, by="gene")
## 3.外显子长度求和
exonic.gene.sizes <- sum(width(reduce(exonic)))
## 4.结果整理
mygeneids <- data.frame(gene_id=names(exonic.gene.sizes),width=exonic.gene.sizes)
