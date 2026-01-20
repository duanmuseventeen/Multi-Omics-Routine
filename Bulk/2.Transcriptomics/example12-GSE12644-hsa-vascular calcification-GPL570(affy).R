rm(list = ls());gc()

# Load pkgs---------------------------------------------------------------------
library(limma)
library(dplyr)
library(stringr)
library(GEOquery)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(affy)

mytheme <- theme_bw() +
  theme(text = element_text(size = 12), 
        plot.margin = ggplot2::margin(30,30,30,30), 
        panel.grid = element_blank(), 
        legend.title = element_blank(),
        legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = "right")

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
# GSE12644 2009 GPL570 ---------------------------------------------------------
# load data --------------------------------------------------------------------
# data are downloaned manulaly from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE12644&format=file
eset <- ReadAffy(celfile.path = "GSE12644_RAW/") %>% 
  rma()
# mas5()

expr_GSE12644 <- eset %>% exprs
# annotate--------------------------------------------------------------------
GPL570 <- "GPL570" %>%
  getGEO()
GPL570_anno <- Table(GPL570)

expr_GSE12644 <- expr_GSE12644 %>% 
  as.data.frame %>%
  mutate(mean = rowMeans(expr_GSE12644)) %>%
  mutate(ID = rownames(expr_GSE12644)) %>%
  left_join(GPL570_anno %>% dplyr::select(ID ,`Gene Symbol`, `ENTREZ_GENE_ID`), by = "ID")  %>%
  filter(`Gene Symbol` != "")
dim(expr_GSE12644)
# [1] 45782    24

expr_GSE12644 <- expr_GSE12644 %>% 
  filter(!str_detect(`Gene Symbol`, "///")) %>% 
  filter(`Gene Symbol` != "---")
dim(expr_GSE12644)
# [1] 42986    24

tmp <- data.frame(
  SYMBOL = unique(expr_GSE12644$`Gene Symbol`),
  stringsAsFactors = F, check.names = F) %>%
  arrange(SYMBOL) %>%
  filter(SYMBOL != "") %>%
  filter(str_sub(SYMBOL,1,3) != "---")

expr_GSE12644_final <- apply(tmp, 1, function(x){
  expr_GSE12644 %>% 
    filter(`Gene Symbol` == x) %>% 
    dplyr::select(contains('CEL')) %>% 
    colMeans(na.rm = T)
})
colnames(expr_GSE12644_final) <- tmp$SYMBOL

head(expr_GSE12644_final[,1:6])
                  A1BG A1BG-AS1     A1CF      A2M  A2M-AS1    A2ML1
GSM317342.CEL 5.274730 5.715667 4.058786 8.996107 6.955682 3.406463
GSM317343.CEL 5.184811 5.892503 4.237803 8.001293 6.872179 3.311361
GSM317344.CEL 5.378121 5.791846 4.112189 7.844297 6.759390 3.305560
GSM317345.CEL 5.536913 5.357080 3.911266 8.135298 7.425385 3.614676
GSM317346.CEL 5.361668 5.856058 3.991763 8.112251 7.038691 3.440475
GSM317347.CEL 5.477980 5.985942 4.174361 8.117654 6.674070 3.309321
# load meta --------------------------------------------------------------------
meta_GSE12644 <- readxl::read_excel("GSE12644_RAW/GSE12644_series_matrix.xlsx")

rownames(expr_GSE12644_final) <- expr_GSE12644_final %>% 
  rownames %>% 
  str_extract("^GSM[0-9]*")

all(rownames(expr_GSE12644_final) == meta_GSE12644$geo_accession)

save(meta_GSE12644, GPL570_anno, expr_GSE12644, expr_GSE12644_final,
     file = "GSE12644_RAW/GSE12644_GPL570_A.Rdata")
# GSVA -------------------------------------------------------------------------
require(GSVA)
require(dplyr)

data <- expr_GSE12644_final %>% t %>% as.data.frame
dim(data)
# 24442    20

hallmarks <-
  c("GOBP_GLYCOSYLATION.v2025.1.Hs.gmt",                 
    "GOBP_PROTEIN_DEGLYCOSYLATION.v2025.1.Hs.gmt",       
    "GOBP_PROTEIN_N_LINKED_GLYCOSYLATION.v2025.1.Hs.gmt",
    "GOBP_PROTEIN_O_LINKED_GLYCOSYLATION.v2025.1.Hs.gmt")
res <- list()
n <- 1
for (i in hallmarks) {
  hallmark <- qusage::read.gmt(i) #返回的是表格
  # old version
  # gsva_matrix <- GSVA::gsva(as.matrix(count), hallmarks, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  
  # new version
  gsvaPar <- ssgseaParam(exprData = as.matrix(data), 
                         geneSets = hallmark,
                         normalize = TRUE)
  gsva_data <- gsva(gsvaPar, verbose = FALSE)
  res[[n]] <- gsva_data
  n <- n + 1
}

res.df <- rbind(res[[1]], res[[2]]) %>% 
  rbind(res[[3]]) %>% 
  rbind(res[[4]])

all(colnames(res.df) == meta_GSE12644$geo_accession)
# visualization-----------------------------------------------------------------
annotation_cols <- data.frame(
  Group = meta_GSE12644$group %>% 
    factor(levels = c("Normal aortic valve","Calcified aortic valve")),
  row.names = colnames(gsva_data)
)
pheatmap::pheatmap(res.df, 
                   clustering_method = "ward.D",
                   annotation_col = annotation_cols,
                   annotation_colors = list(
                     Group = c(`Calcified aortic valve` = "#C00000", 
                               `Normal aortic valve` = "black")),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "row")
  


