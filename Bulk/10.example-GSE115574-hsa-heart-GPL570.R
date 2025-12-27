rm(list = ls());gc()

# Load pkgs---------------------------------------------------------------------
library(limma)
library(stringr)
library(GEOquery)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(affy)
# my functions------------------------------------------------------------------
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

# GSE115574---------------------------------------------------------------------
if(GSE115574){
  
  # because of poor connection, data are downloaned manulaly from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115574&format=file
  eset <- ReadAffy(celfile.path = "GSE115574/GSE115574_RAW/") %>% 
    rma()
  # mas5()

  # some plaform cannot be read by affy pkg, pls use oligo pkg according to warning
  # library(oligo)
  # eset <- oligo::read.celfiles(filenames = dir("GSE77287_RAW/", full.names = T)) %>% 
  #  oligo::rma()
  
  expr_GSE115574 <- eset %>% exprs
  # annotate--------------------------------------------------------------------
  GPL570 <- "GPL570" %>%
    getGEO()
  GPL570_anno <- Table(GPL570)
  
  expr_GSE115574 <- expr_GSE115574 %>% 
    as.data.frame %>%
    mutate(mean = rowMeans(expr_GSE115574)) %>%
    mutate(ID = rownames(expr_GSE115574)) %>%
    left_join(GPL570_anno %>% dplyr::select(ID ,`Gene Symbol`, `ENTREZ_GENE_ID`), by = "ID")  %>%
    filter(`Gene Symbol` != "")
  dim(expr_GSE115574)
  [1] 45782    63
  
  for (i in 1:nrow(expr_GSE115574)) {
    tmp <- str_split(expr_GSE115574$`Gene Symbol`[i],"///")
    if(length(tmp[[1]]) > 1){
      for (k in 1:length(tmp[[1]])) {
        tmp[[1]][k] <- tmp[[1]][k] %>% str_remove_all(" ")
      }
      tmp[[1]] <- tmp[[1]][which(tmp[[1]] != "---")]
      if(length(tmp[[1]]) == 1){
        expr_GSE115574$`Gene Symbol`[i] <- tmp[[1]][1]
      }else if(length(tmp[[1]]) > 1){
        tmprow <- expr_GSE115574[i,]
        for (j in 1:(length(tmp[[1]]) - 1)) {
          expr_GSE115574 <- rbind(expr_GSE115574, tmprow)
          expr_GSE115574$`Gene Symbol`[nrow(expr_GSE115574)] <- tmp[[1]][[j + 1]]
        }
        expr_GSE115574$`Gene Symbol`[i] <- tmp[[1]][[1]]
      }
    }
  }
  
  dim(expr_GSE115574)
  [1] 50362    63
  tmp <- data.frame(
    SYMBOL = unique(expr_GSE115574$`Gene Symbol`),
    stringsAsFactors = F, check.names = F) %>%
    arrange(SYMBOL) %>%
    filter(SYMBOL != "") %>%
    filter(str_sub(SYMBOL,1,3) != "---")
  
  expr_GSE115574_final <- apply(tmp, 1, function(x){
    expr_GSE115574 %>% filter(`Gene Symbol` == x) %>% 
      dplyr::select(eset %>% exprs %>% colnames) %>% colMeans(na.rm = T)
    })
  colnames(expr_GSE115574_final) <- tmp$SYMBOL
  
  meta_GSE115574 <- readxl::read_excel("GSE115574/GSE115574_series_matrix.xlsx")
  
  rownames(expr_GSE115574_final) <- expr_GSE115574_final %>% 
    rownames %>% 
    str_extract("^GSM[0-9]*")
  
  all(rownames(expr_GSE115574_final) == meta_GSE115574$geo_accession)
  
  save(meta_GSE115574, GPL570_anno, expr_GSE115574, expr_GSE115574_final,
       file = "GSE115574/GSE115574_GPL570_A.Rdata")
  
  # Exploration-----------------------------------------------------------------
  search <- c()
  # all
  res.all <- purrr::map(search, function(x){
    t.test(expr_GSE115574_final[,x] ~ meta_GSE115574$disease) 
  })

  # left atrium
  meta.left <- meta_GSE115574 %>% filter(tissue == "left atrium  - heart")
  expr.left <- expr_GSE115574_final[meta.left$geo_accession,]
  res.left <- purrr::map(search, function(x){
    t.test(expr.left[,x] ~ meta.left$disease) 
  })
  
  # right atrium
  meta.right <- meta_GSE115574 %>% filter(tissue == "right atrium  - heart")
  expr.right <- expr_GSE115574_final[meta.right$geo_accession,]
  res.right <- purrr::map(search, function(x){
    t.test(expr.right[,x] ~ meta.right$disease) 
  })
}

