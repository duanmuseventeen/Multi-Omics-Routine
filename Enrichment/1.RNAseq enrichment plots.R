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
library(DEXSeq)

mytheme <- theme_bw() +
  theme(text = element_text(size = 12), 
        plot.margin = ggplot2::margin(30,30,30,30), 
        panel.grid = element_blank(), 
        legend.title = element_blank(),
        legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = "right")

# ERA enrichment----------------------------------------------------------------
# In general, geneList are names of DEGs 
ERA_KEGG <- enrichKEGG(
  gene = geneList,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1, 
  qvalueCutoff = 1)
ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")

# for Mus musculus only
ERA_KEGG@result$Description <- ERA_KEGG@result$Description %>% 
  str_remove(" - Mus musculus \\(house mouse\\)") %>% 
  dotplot(showCategory = 20) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

ERA_GO <- enrichGO(
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
ERA_GO <- setReadable(ERA_GO, 'org.Mm.eg.db', 'ENTREZID')

BPtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "BP");BPtop10 <- BPtop10$ID[1:10]
CCtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "CC");CCtop10 <- CCtop10$ID[1:10]
MFtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "MF");MFtop10 <- MFtop10$ID[1:10]

ERA_GO@result %>% 
  filter(ID %in% c(BPtop10, CCtop10, MFtop10)) %>% 
  dotplot(
    col = "pvalue", 
    showCategory = 50) + 
    facet_grid(ONTOLOGY~., scale='free') + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))

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

# geneList2 are log2FC named with entrezID
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