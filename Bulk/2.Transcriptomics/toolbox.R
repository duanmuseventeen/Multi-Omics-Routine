# main -------------------------------------------------------------------------
suppressWarnings(require(crayon, warn.conflicts = F))
suppressWarnings(require(dplyr, warn.conflicts = F))
suppressWarnings(require(stringr, warn.conflicts = F))
suppressWarnings(require(DESeq2, warn.conflicts = F))
suppressWarnings(require(clusterProfiler, warn.conflicts = F))
suppressWarnings(require(biomaRt, warn.conflicts = F))
suppressMessages(require(ggplot2))
suppressMessages(require(FactoMineR))
suppressMessages(require(limma))
suppressMessages(require(sva))

Pa    = crayon::cyan
Er    = crayon::red$bold
Sa    = crayon::blue
No    = crayon::magenta$bold
Wa    = crayon::yellow

checkPath = function(x) normalizePath(x, '/', T)


myORA <- function(data = DEG){
  t0 <- Sys.time()
  cat("\n[", format(t0, "%Y-%m-%d %H:%M:%S"), "] ORA START\n", sep = "")
  
  set.seed(1011)
  
  DEG4enrichment <- data %>% 
    mutate(hgnc_symbol = rownames(data)) %>% 
    left_join(conversion_table, by = "hgnc_symbol") %>% 
    mutate(sig = case_when(abs(log2FoldChange) >= 1 & padj < .05 ~ 1,
                           TRUE ~ 0)) %>% 
    filter(sig == 1)
  
  ORA_KEGG <- enrichKEGG(
    gene = DEG4enrichment$entrezgene_id,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1, 
    qvalueCutoff = 1)
  ORA_KEGG <- setReadable(ORA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  ORA_GO <- enrichGO(
    gene = DEG4enrichment$entrezgene_id,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID", 
    ont = "ALL",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = T)
  ORA_GO <- setReadable(ORA_GO, 'org.Hs.eg.db', 'ENTREZID')
  
  ORA_Reactome <- ReactomePA::enrichPathway(
    gene= DEG4enrichment$entrezgene_id, 
    pvalueCutoff = 1, 
    readable=TRUE, 
    organism = "human")
  
  t1 <- Sys.time()
  cat("[", format(t1, "%Y-%m-%d %H:%M:%S"), "] ORA END; elapsed = ",
      round(as.numeric(difftime(t1, t0, units = "mins")), 2), " mins\n", sep = "")
  
  ORA = list(KEGG = ORA_KEGG, GO = ORA_GO, REACTOME = ORA_Reactome)
}

mygsea <- function(data = DEG, nPerm = 100000){
  t0 <- Sys.time()
  cat("\n[", format(t0, "%Y-%m-%d %H:%M:%S"), "] GSEA START\n", sep = "")
  
  set.seed(1011)
  
  DEG4enrichment <- data %>% 
    mutate(hgnc_symbol = rownames(data)) %>% 
    left_join(conversion_table, by = "hgnc_symbol") %>% 
    filter(complete.cases(entrezgene_id)) %>% 
    arrange(desc(abs(log2FoldChange))) %>% 
    filter(!duplicated(entrezgene_id)) %>% 
    arrange(desc(log2FoldChange))
  
  geneList2 <- DEG4enrichment$log2FoldChange
  names(geneList2) <- DEG4enrichment$entrezgene_id
  geneList2 <- geneList2 + sort(rnorm(length(geneList2), mean = 0, sd = 1e-9), decreasing = TRUE)
  
  GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                       organism = "hsa",
                       minGSSize    = 10,
                       pvalueCutoff = 1,
                       verbose      = FALSE, 
                       nPermSimple = nPerm,
                       eps = 0)
  warnings() %>% print
  GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  GSEA_GO <- gseGO(geneList = geneList2,
                   ont = "ALL",
                   OrgDb = "org.Hs.eg.db",
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = FALSE, 
                   nPermSimple = nPerm,
                   eps = 0)
  warnings() %>% print
  GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  GSEA_REACTOME <- ReactomePA::gsePathway(geneList2,
                                          organism= "human",
                                          minGSSize= 10,
                                          maxGSSize= 500,
                                          pvalueCutoff= 1,
                                          pAdjustMethod= "BH",
                                          verbose= FALSE,
                                          nPermSimple = nPerm,
                                          eps= 0)
  warnings() %>% print
  GSEA_REACTOME <- setReadable(GSEA_REACTOME, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  t1 <- Sys.time()
  cat("[", format(t1, "%Y-%m-%d %H:%M:%S"), "] GSEA END; elapsed = ",
      round(as.numeric(difftime(t1, t0, units = "mins")), 2), " mins\n", sep = "")
  
  GSEA <- list(KEGG = GSEA_KEGG, GO = GSEA_GO, REACTOME = GSEA_REACTOME)
}

myhallmark1 <- function(data = DEG, nPerm = 100000){
  t0 <- Sys.time()
  cat("\n[", format(t0, "%Y-%m-%d %H:%M:%S"), "] GSEA HALLMARK START\n", sep = "")
  
  geneset <- read.gmt("E:/cold2hot/h.all.v2026.1.Hs.symbols.gmt")
  
  geneList2 <- data$log2FoldChange
  names(geneList2) <- rownames(data)
  geneList2 <- geneList2 + sort(rnorm(length(geneList2), mean = 0, sd = 1e-9), decreasing = TRUE)
  
  set.seed(1011)
  gsea_h50 <- GSEA(geneList2, 
                   TERM2GENE = geneset,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   minGSSize = 10,  
                   maxGSSize = 500, 
                   nPermSimple = nPerm,
                   verbose = FALSE)
  
  t1 <- Sys.time()
  cat("[", format(t1, "%Y-%m-%d %H:%M:%S"), "] GSEA HALLMARK END; elapsed = ",
      round(as.numeric(difftime(t1, t0, units = "mins")), 2), " mins\n", sep = "")
  
  return(gsea_h50)
}

mysave <- function(datasetid){
  save(DEG,  file = paste0(datasetid,'-DEGs.Rdata'))
  if(sum(abs(DEG$log2FoldChange) >= 1 & DEG$padj < .05) > 10) save(ORA, file = paste0(datasetid,'-ORA.Rdata'))
  save(GSEA, file = paste0(datasetid,'-GSEA.Rdata'))
  write.csv(DEG, paste0(datasetid,"-DEG.csv"))
  if(sum(abs(DEG$log2FoldChange) >= 1 & DEG$padj < .05) > 10)  for(i in 1:3){write.csv(ORA [[i]]@result, paste0(datasetid, "-ORA_", names(ORA)[i],".csv"))}
  for(i in 1:3){write.csv(GSEA[[i]]@result, paste0(datasetid, "-GSEA_",names(GSEA)[i],".csv"))}
  save(gsea_h50, file = paste0(datasetid,"-GSEA_HALLMARK50.Rdata"))
}

mypca <- function(data, features, mytitle, plot = "group"){
  features <- intersect(features, rownames(data))
  pre.pca <- PCA(data[features,] %>% t, graph = FALSE)
  if(plot == "group"){
    factoextra::fviz_pca_ind(pre.pca,
                             mean.point = F,
                             pointsize = 1,
                             geom= "point",
                             col.ind = meta$group,
                             addEllipses = T) +
      guides(size = "none", shape = "none", fill = "none") +
      scale_color_manual(values= c("#4798ac", "#b5333b")) +
      scale_fill_manual(values = c("#4798ac", "#b5333b")) +
      labs(col =  "Cancer type", title = mytitle) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 12),
            plot.margin =  ggplot2::margin(12,12,12,12), panel.grid = element_blank(),
            legend.background = element_rect(linetype = 1, color = "black"), legend.position = "right",
            axis.text = element_text(size= 16, color = "black"))
  }else if(plot == "purity"){
    factoextra::fviz_pca_ind(pre.pca,
                             mean.point = F,
                             pointsize = 1,
                             geom= "point",
                             col.ind = meta$purity,
                             addEllipses = F) +
      guides(size = "none", shape = "none", fill = "none") +
      scale_color_gradient2(low = "#4DBBD5FF", mid = "gray80", high = "#E64B35FF", midpoint = 0.5) +
      labs(col =  "Cancer type", title = mytitle) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 12),
            plot.margin =  ggplot2::margin(12,12,12,12), panel.grid = element_blank(),
            legend.background = element_rect(linetype = 1, color = "black"), legend.position = "right",
            axis.text = element_text(size= 16, color = "black"))
  }
}