# Reference
# https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/
# https://guangchuangyu.github.io/2016/07/leading-edge-analysis/
# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
# https://github.com/kerseviciute/aPEAR
# https://mp.weixin.qq.com/s/eX8mc5ssTQht8rhU0bbp6g
# https://www.jianshu.com/p/de8fccabc2e7
# https://zhuanlan.zhihu.com/p/830574147

library(stringr)
library(dplyr)
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(enrichplot)
require(survminer)
require(survival)
library(clusterProfiler)
library(ReactomePA)
# library(GseaVis)
# devtools::install_github("nicolash2/gggsea")
# library(gggsea)

mytheme <- theme_bw() +
  theme(text = element_text(size = 12), 
        plot.margin = ggplot2::margin(30,30,30,30), 
        panel.grid = element_blank(), 
        legend.title = element_blank(),
        legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = "right")

# ORA enrichment----------------------------------------------------------------
# In general, geneList are names of DEGs 
ORA_KEGG <- enrichKEGG(
  gene = geneList,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1, 
  qvalueCutoff = 1)
ORA_KEGG <- setReadable(ORA_KEGG, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")

# for Mus musculus only
ORA_KEGG@result$Description <- ORA_KEGG@result$Description %>% 
  str_remove(" - Mus musculus \\(house mouse\\)")
dotplot(ORA_KEGG, showCategory = 20) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

ORA_GO <- enrichGO(
  gene = geneList,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENTREZID", 
  ont = "ALL",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = T)
ORA_GO <- setReadable(ORA_GO, 'org.Mm.eg.db', 'ENTREZID')

BPtop10 <- ORA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "BP");BPtop10 <- BPtop10$ID[1:10]
CCtop10 <- ORA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "CC");CCtop10 <- CCtop10$ID[1:10]
MFtop10 <- ORA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "MF");MFtop10 <- MFtop10$ID[1:10]

ORA_GO@result %>% 
  filter(ID %in% c(BPtop10, CCtop10, MFtop10))
                   
dotplot(
    ORA_GO,
    col = "pvalue", 
    showCategory = 50) + 
    facet_grid(ONTOLOGY~., scale='free') + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))

ORA_Reactome <- ReactomePA::enrichPathway(
  gene= geneList, 
  pvalueCutoff = 1, 
  readable=TRUE, 
  organism = "human")
clusterProfiler::dotplot(ORA_Reactome)

# GSEA enrichment---------------------------------------------------------------
# NES can be DEA results from gseKEGG or gseGO
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

gg_KEGGE_dot <- function(NES, Topn = 10){
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
  ggplot(tmp, aes(x = NES, y = Description, color = NES, size = -log10(p.adjust))) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_point() +
    scale_color_gradient2(low = "#13679E", mid = "#c4a755", high = "#AB3A29") +
    scale_size_continuous(range=c(2,6)) +
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

gg_GOE_dot <- function(NES, Topn = 10){
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
  ggplot(tmp, aes(x = NES, y = ID, color = NES, size = -log10(p.adjust))) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_point() +
    scale_color_gradient2(low = "#13679E", mid = "#c4a755", high = "#AB3A29") +
    scale_size_continuous(range=c(2,6)) +
    labs(y = "") +
    guides() +
    facet_grid(ONTOLOGY~., scale='free') + 
    mytheme
}

# geneList2 are sorted log2FC named with entrezID
GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                     seed = 1011,
                     organism = "mmu",
                     minGSSize    = 10,
                     pvalueCutoff = 1,
                     verbose      = FALSE, 
                     eps = 0)
GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")

GSEA_GO <- gseGO(geneList = genes,
                 ont = "ALL",
                 OrgDb = "org.Mm.eg.db",
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 verbose      = FALSE, 
                 eps = 0)
GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")

# GSEA with fgsea-------------------------------------------------------------                   
# geneList2 are sorted log2FC with names
library(fgsea)

gmt <- qusage::read.gmt("HP_HEPATIC_FIBROSIS.v2024.1.Hs.gmt") 
fgseaRes <- fgsea(gmt, geneList2, minSize = 5, maxSize = 500, nperm=1000, gseaParam = 1)
# pathway      pval      padj    log2err        ES       NES  size
# <char>     <num>     <num>      <num>     <num>     <num> <int>
# 1: HP_HEPATIC_FIBROSIS 0.9775785 0.9775785 0.05195125 0.1966983 0.7319902   122
#    leadingEdge
#         <list>
# 1: ARHGAP31....

plotEnrichment(gmt$HP_HEPATIC_FIBROSIS, geneList2)
plotGseaTable(gmt, geneList2, fgseaRes, gseaParam=0.5)
# Visualization--------------------------------------------------------------- 
require(aPEAR)
p <- enrichmentNetwork(as.data.frame(Reactome@result), drawEllipses = TRUE,
                  fontSize = 4, repelLabels = TRUE)

library(plotly)
ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
