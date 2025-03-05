library(limma)
library(string)
library(GEOquery)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

mytheme <- theme_bw() +
  theme(text = element_text(size = 12), 
        plot.margin = ggplot2::margin(30,30,30,30), 
        panel.grid = element_blank(), 
        legend.title = element_blank(),
        legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = "right")

pca <- function(dat, group){
  
  suppressMessages(require(ggplot2))
  suppressMessages(require(FactoMineR))
  
  pre.pca <- PCA(dat, graph = FALSE)
  factoextra::fviz_pca_ind(pre.pca,
                           mean.point = F,
                           pointsize = 4,
                           geom= "point",
                           col.ind = group,
                           addEllipses = F) +
    stat_ellipse(level = 0.95, type = "norm", geom = "polygon", col = "gray20", alpha = 0) +
    # geom_point(size = 4) +
    # ggtitle("") +
    # scale_shape_manual(values = ) +
    # scale_fill_manual(values = ) +
    # scale_color_manual(values = ) +
    guides(size = "none") +
    labs(col = "Group", shape = "Group") +
    ggtitle("") +
    theme_bw() +
    theme(text = element_text(size = 20),
          plot.margin =  ggplot2::margin(30,30,30,30),
          panel.grid = element_blank(),
          legend.background = element_rect(linetype = 1, color = "black"),
          legend.position = "right",
          axis.text = element_text(size= 16, color = "black"))
}

gg_KEGGE_bar <- function(NES, Topn = 10){
  # for Mus musculus only
  NES@result$Description <- NES@result$Description %>% str_remove(" - Mus musculus \\(house mouse\\)") 
  
  Top <- NES@result
  TopUp <- Top %>% 
    filter(p.adjust < .05 & NES > 0) %>%
    arrange(desc(NES))
  
  TopDown <- Top %>% 
    filter(p.adjust < .05 & NES < 0) %>%
    arrange(NES)
  
  KEGG_EA2 <- NES
  KEGG_EA2@result <- KEGG_EA2@result %>% 
    filter(Description %in% c(TopUp$Description[1:Topn], TopDown$Description[1:Topn]))
  
  tmp <- KEGG_EA2@result %>% arrange(NES) %>%
    mutate(Description = factor(Description, levels = Description),
           color = ifelse(NES < 0, 0, 1) %>% factor(labels = c("down","up")))
  ggplot(tmp, aes(x = NES, y = Description, fill = color)) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("#448844","#AB3A29")) +
    labs(y = "") +
    guides(fill = "none") +
    mytheme
}

gg_GOE_bar <- function(NES, Topn = 10){
  Top <- NES@result
  
  TopBPUp <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "BP" & NES > 0) %>%
    arrange(desc(NES))
  
  TopBPDown <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "BP" & NES < 0) %>%
    arrange(NES)
  
  TopCCUp <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "CC" & NES > 0) %>%
    arrange(desc(NES))
  
  TopCCDown <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "CC" & NES < 0) %>%
    arrange(NES)
  
  TopMFUp <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "MF" & NES > 0) %>%
    arrange(desc(NES))
  
  TopMFDown <- Top %>% 
    filter(p.adjust < .05 & ONTOLOGY == "MF" & NES < 0) %>%
    arrange(NES)
  
  GO_EA2 <- NES
  GO_EA2@result <- GO_EA2@result %>% 
    filter(Description %in% c(TopBPUp$Description[1:Topn], TopBPDown$Description[1:Topn],
                              TopCCUp$Description[1:Topn], TopCCDown$Description[1:Topn],
                              TopMFUp$Description[1:Topn], TopMFDown$Description[1:Topn]))
  
  tmp <- GO_EA2@result %>% arrange(NES) %>%
    mutate(ID = factor(ID, levels = ID),
           color = ifelse(NES < 0, 0, 1) %>% factor(labels = c("down","up")))
  ggplot(tmp, aes(x = NES, y = ID, fill = color)) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("#448844","#AB3A29")) +
    labs(y = "") +
    guides(fill = "none") +
    facet_grid(ONTOLOGY~., scale='free') + 
    mytheme
}

# GSE84044----------------------------------------------------------------------
if(GSE84044){
  
  # get data--------------------------------------------------------------------
  GSE84044 <- getGEO("GSE84044")
  # GSE84044 %>% length
  
  meta_GSE84044 <- GSE84044[[1]] %>% pData()
  expr_GSE84044 <- GSE84044[[1]] %>% exprs
  dim(expr_GSE84044)
  [1] 54675   124
  
  # annotate--------------------------------------------------------------------
  GPL570 <- GSE84044[[1]]$platform_id[1] %>% 
    as.character() %>%
    getGEO()
  GPL570_anno <- Table(GPL570)
  
  expr_GSE84044 <- expr_GSE84044 %>% 
    as.data.frame %>%
    mutate(mean = rowMeans(expr_GSE84044)) %>%
    mutate(ID = rownames(expr_GSE84044)) %>%
    left_join(GPL570_anno %>% dplyr::select(ID ,`Gene Symbol`, `ENTREZ_GENE_ID`), by = "ID")  %>%
    filter(`Gene Symbol` != "")
  dim(expr_GSE84044)
  [1] 45782   128
  
  for (i in 1:nrow(expr_GSE84044)) {
    tmp <- str_split(expr_GSE84044$`Gene Symbol`[i],"///")
    if(length(tmp[[1]]) > 1){
      for (k in 1:length(tmp[[1]])) {
        tmp[[1]][k] <- tmp[[1]][k] %>% str_remove_all(" ")
      }
      tmp[[1]] <- tmp[[1]][which(tmp[[1]] != "---")]
      if(length(tmp[[1]]) == 1){
        expr_GSE84044$`Gene Symbol`[i] <- tmp[[1]][1]
      }else if(length(tmp[[1]]) > 1){
        tmprow <- expr_GSE84044[i,]
        for (j in 1:(length(tmp[[1]]) - 1)) {
          expr_GSE84044 <- rbind(expr_GSE84044, tmprow)
          expr_GSE84044$`Gene Symbol`[nrow(expr_GSE84044)] <- tmp[[1]][[j + 1]]
        }
        expr_GSE84044$`Gene Symbol`[i] <- tmp[[1]][[1]]
      }
    }
  }
  dim(expr_GSE84044)
  [1] 50362    76
  
  tmp <- data.frame(
    SYMBOL = unique(expr_GSE84044$`Gene Symbol`),
    stringsAsFactors = F, check.names = F) %>%
    arrange(SYMBOL) %>%
    filter(SYMBOL != "") %>%
    filter(str_sub(SYMBOL,1,3) != "---")
  expr_GSE84044_final <- apply(tmp, 1, function(x){
    expr_GSE84044 %>% filter(`Gene Symbol` == x) %>% 
      dplyr::select(colnames(exprs(GSE84044[[1]]))) %>% colMeans(na.rm = T)})
  colnames(expr_GSE84044_final) <- tmp$SYMBOL
  save(GSE84044, meta_GSE84044, GPL570_anno, expr_GSE84044, expr_GSE84044_final,
       file = "GSE84044/GSE84044_GPL570_A.Rdata")
  
  # DEA by limma----------------------------------------------------------------
  # Because the Experiment type	of GSE84044 is 'Expression profiling by array',
  # the limma pkg is selected to perform DEA.
  
  # Liver biopsy is considered to be the “gold standard” for the diagnosis of fibrosis, and sequential histological staging of fibrosis (Scheuer score ‘S’) and grading of inflammation (Scheuer score ‘G’) can be used to assess disease progression.
  
  array <- expr_GSE84044_final %>% t %>% as.data.frame
  group_list <- meta_GSE84044$`scheuer score s:ch1`
  batch <- meta_GSE84044$`batch:ch1`
  
  pca(t(array), batch) # no apparent batch effect
  pca(t(array), group_list) # the difference among 5 groups are not signifcant
  
  names(group_list) <- meta_GSE84044$geo_accession
  group_list <- group_list %>% as.numeric() 
  group_list[group_list < 2] <- 0
  group_list[group_list >=2] <- 1
  group_list <- group_list %>% factor(labels = c("Normal", "Fibrosis"))
  
  all(names(group_list) == colnames(array))
  # [1] TRUE
  
  design <- model.matrix(~0+group_list)
  rownames(design) = colnames(array)
  colnames(design) <- levels(group_list)
  cont.matrix <- makeContrasts(contrasts = c('Fibrosis-Normal'), levels = design)
  
  fit <- lmFit(array, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  degs <-topTable(fit2, coef = 'Fibrosis-Normal', n = Inf)
  
  # Enrichment------------------------------------------------------------------
  # ORA
  genelist <- bitr(rownames(degs)[degs$adj.P.Val < 0.05 & abs(degs$logFC) >= 0.5849625],
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  
  ERA_KEGG <- enrichKEGG(
    gene = genelist$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1, 
    qvalueCutoff = 1)
  ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  dotplot(ERA_KEGG, showCategory = 40)
  
  ERA_GO <- enrichGO(
    gene = genelist$ENTREZID,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID", 
    ont = "ALL",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = T)
  ERA_GO <- setReadable(ERA_GO, 'org.Hs.eg.db', 'ENTREZID')
  
  BPtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "BP");BPtop10 <- BPtop10$ID[1:10]
  CCtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "CC");CCtop10 <- CCtop10$ID[1:10]
  MFtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "MF");MFtop10 <- MFtop10$ID[1:10]
  
  ERA_GO_R <- ERA_GO
  ERA_GO_R@result <- ERA_GO_R@result %>% 
    filter(ID %in% c(BPtop10, CCtop10, MFtop10))
  dotplot(ERA_GO_R, col = "pvalue", showCategory = 50) + 
    facet_grid(ONTOLOGY~., scale='free') + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  
  # GSEA
  degs  <- degs %>% arrange(desc(logFC))
  annot <- bitr(rownames(degs), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  annot <- annot %>% 
    left_join(degs %>% tibble::rownames_to_column("SYMBOL"), by = "SYMBOL")
  genes <- annot$logFC
  names(genes) <- annot$ENTREZID
  
  GSEA_KEGG <- gseKEGG(geneList = genes, 
                       seed = 1011,
                       organism = "hsa",
                       minGSSize    = 10,
                       pvalueCutoff = 1,
                       verbose      = FALSE, 
                       eps = 0)
  GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  gg_KEGGE_bar(GSEA_KEGG, Topn = 30)
  
  which(GSEA_KEGG@result$ID %in% c('hsa04630')) # JAK-STAT
  p_KEGG_GSEA <- enrichplot::gseaplot2(GSEA_KEGG, c(120), pvalue_table = TRUE)
  
  GSEA_GO <- gseGO(geneList = genes,
                   ont = "ALL",
                   OrgDb = "org.Hs.eg.db",
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = FALSE, 
                   eps = 0)
  GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  gg_GOE_bar(GSEA_GO)
  
  # analysis of cytokines-------------------------------------------------------------------
  # the list of totally 84 cytokine genes were retrievaled from 38057668
  cytokines <- readxl::read_excel("cytokines from 38057668.xlsx")
  cytokines <- cytokines %>% 
    filter(complete.cases(Symbol))
  cytokine <- cytokines$Symbol %>% unique
  
  array.cytokines <- array %>% 
    tibble::rownames_to_column("Symbol") %>% 
    filter(Symbol %in% cytokine) %>% 
    tibble::column_to_rownames("Symbol") %>% 
    as.data.frame
  dim(array.cytokines)
  [1]  82 124
  
  names(group_list) <- meta_GSE84044$geo_accession
  group_list <- group_list %>% factor(labels = c("Normal", "Fibrosis"))
  
  all(names(group_list) == colnames(array))
  
  annotation_col = data.frame(
    `Scheuer Score (S)` = meta_GSE84044$`scheuer score s:ch1`,
    row.names = meta_GSE84044$geo_accession,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>% 
    mutate(group = ifelse(`Scheuer Score (S)` %in% c("0", "1"), 0, 1) %>% 
             factor(labels = c("Normal", "Fibrosis")))
  
  annotation_col <- annotation_col %>% arrange(`Scheuer Score (S)`) 
  array.cytokines <- array.cytokines[rownames(annotation_col)]
  
  all(rownames(annotation_col) == colnames(array.cytokines))
  
  pheatmap::pheatmap(array.cytokines, 
                     color = colorRampPalette(c("white","red3"),alpha = TRUE)(200),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     annotation_col = annotation_col,
                     clustering_method = "ward.D2",
                     scale = "row")
}
