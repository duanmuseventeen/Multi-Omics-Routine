qc.dist.heat <- function(dat){
  
  suppressMessages(require(dplyr))
  suppressMessages(require(RColorBrewer))
  suppressMessages(require(pheatmap))
  
  # dat with rowname should only includes QC samples, row means sample, column means feature
  sampleDists <- dist(dat)
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  rownames(sampleDistMatrix) <- rownames(dat)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col=colors)
}

