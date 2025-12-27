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

# GSE79768---------------------------------------------------------------------
if(GSE79768){
  
  # because of poor connection, data are downloaned manulaly from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE79768&format=file
  eset <- ReadAffy(celfile.path = "GSE79768/GSE79768_RAW/") %>% 
    rma()
  # mas5()
  
  expr_GSE79768 <- eset %>% exprs
  # annotate--------------------------------------------------------------------
  GPL570 <- "GPL570" %>%
    getGEO()
  GPL570_anno <- Table(GPL570)
  
  expr_GSE79768 <- expr_GSE79768 %>% 
    as.data.frame %>%
    mutate(mean = rowMeans(expr_GSE79768)) %>%
    mutate(ID = rownames(expr_GSE79768)) %>%
    left_join(GPL570_anno %>% dplyr::select(ID ,`Gene Symbol`, `ENTREZ_GENE_ID`), by = "ID")  %>%
    filter(`Gene Symbol` != "")
  dim(expr_GSE79768)
  [1] 45782    30
  
  for (i in 1:nrow(expr_GSE79768)) {
    tmp <- str_split(expr_GSE79768$`Gene Symbol`[i],"///")
    if(length(tmp[[1]]) > 1){
      for (k in 1:length(tmp[[1]])) {
        tmp[[1]][k] <- tmp[[1]][k] %>% str_remove_all(" ")
      }
      tmp[[1]] <- tmp[[1]][which(tmp[[1]] != "---")]
      if(length(tmp[[1]]) == 1){
        expr_GSE79768$`Gene Symbol`[i] <- tmp[[1]][1]
      }else if(length(tmp[[1]]) > 1){
        tmprow <- expr_GSE79768[i,]
        for (j in 1:(length(tmp[[1]]) - 1)) {
          expr_GSE79768 <- rbind(expr_GSE79768, tmprow)
          expr_GSE79768$`Gene Symbol`[nrow(expr_GSE79768)] <- tmp[[1]][[j + 1]]
        }
        expr_GSE79768$`Gene Symbol`[i] <- tmp[[1]][[1]]
      }
    }
  }
  
  dim(expr_GSE79768)
  [1] 50362    63
  tmp <- data.frame(
    SYMBOL = unique(expr_GSE79768$`Gene Symbol`),
    stringsAsFactors = F, check.names = F) %>%
    arrange(SYMBOL) %>%
    filter(SYMBOL != "") %>%
    filter(str_sub(SYMBOL,1,3) != "---")
  
  expr_GSE79768_final <- apply(tmp, 1, function(x){
    expr_GSE79768 %>% filter(`Gene Symbol` == x) %>% 
      dplyr::select(eset %>% exprs %>% colnames) %>% colMeans(na.rm = T)
  })
  colnames(expr_GSE79768_final) <- tmp$SYMBOL
  
  meta_GSE79768 <- readxl::read_excel("GSE79768/GSE79768_series_matrix.xlsx")
  
  rownames(expr_GSE79768_final) <- expr_GSE79768_final %>% 
    rownames %>% 
    str_extract("^GSM[0-9]*")
  
  all(rownames(expr_GSE79768_final) == meta_GSE79768$geo_accession)
  
  save(meta_GSE79768, GPL570_anno, expr_GSE79768, expr_GSE79768_final,
       file = "GSE79768/GSE79768_GPL570_A.Rdata")
  
  # Exploration-----------------------------------------------------------------
  search <- c()

  res.all <- purrr::map(search, function(x){
    t.test(expr_GSE79768_final[,x] ~ meta_GSE79768$disease) 
  })
}