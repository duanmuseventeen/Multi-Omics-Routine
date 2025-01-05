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

# GSE14975---------------------------------------------------------------------
if(GSE14975){
  
  # because of poor connection, data are downloaned manulaly from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14975&format=file
  eset <- ReadAffy(celfile.path = "GSE14975/GSE14975_RAW/") %>% 
    rma()
  # mas5()
  
  expr_GSE14975 <- eset %>% exprs
  # annotate--------------------------------------------------------------------
  GPL570 <- "GPL570" %>%
    getGEO()
  GPL570_anno <- Table(GPL570)
  
  expr_GSE14975 <- expr_GSE14975 %>% 
    as.data.frame %>%
    mutate(mean = rowMeans(expr_GSE14975)) %>%
    mutate(ID = rownames(expr_GSE14975)) %>%
    left_join(GPL570_anno %>% dplyr::select(ID ,`Gene Symbol`, `ENTREZ_GENE_ID`), by = "ID")  %>%
    filter(`Gene Symbol` != "")
  dim(expr_GSE14975)
  [1] 45782    14
  
    for (i in 1:nrow(expr_GSE14975)) {
      tmp <- str_split(expr_GSE14975$`Gene Symbol`[i],"///")
      if(length(tmp[[1]]) > 1){
        for (k in 1:length(tmp[[1]])) {
          tmp[[1]][k] <- tmp[[1]][k] %>% str_remove_all(" ")
        }
        tmp[[1]] <- tmp[[1]][which(tmp[[1]] != "---")]
        if(length(tmp[[1]]) == 1){
          expr_GSE14975$`Gene Symbol`[i] <- tmp[[1]][1]
        }else if(length(tmp[[1]]) > 1){
          tmprow <- expr_GSE14975[i,]
          for (j in 1:(length(tmp[[1]]) - 1)) {
            expr_GSE14975 <- rbind(expr_GSE14975, tmprow)
            expr_GSE14975$`Gene Symbol`[nrow(expr_GSE14975)] <- tmp[[1]][[j + 1]]
          }
          expr_GSE14975$`Gene Symbol`[i] <- tmp[[1]][[1]]
        }
      }
    }
  
  dim(expr_GSE14975)
  [1] 50362    14
  tmp <- data.frame(
    SYMBOL = unique(expr_GSE14975$`Gene Symbol`),
    stringsAsFactors = F, check.names = F) %>%
    arrange(SYMBOL) %>%
    filter(SYMBOL != "") %>%
    filter(str_sub(SYMBOL,1,3) != "---")
  
  expr_GSE14975_final <- apply(tmp, 1, function(x){
    expr_GSE14975 %>% filter(`Gene Symbol` == x) %>% 
      dplyr::select(eset %>% exprs %>% colnames) %>% colMeans(na.rm = T)
  })
  colnames(expr_GSE14975_final) <- tmp$SYMBOL
  
  meta_GSE14975 <- readxl::read_excel("GSE14975/GSE14975_series_matrix.xlsx")
  
  rownames(expr_GSE14975_final) <- expr_GSE14975_final %>% 
    rownames %>% 
    str_extract("^GSM[0-9]*")
  
  all(rownames(expr_GSE14975_final) == meta_GSE14975$geo_accession)
  
  save(meta_GSE14975, GPL570_anno, expr_GSE14975, expr_GSE14975_final,
       file = "GSE14975/GSE14975_GPL570_A.Rdata")
  
  # Exploration-----------------------------------------------------------------
  search <- c("SMPD1","SMPD2","SMPD3")
  
  table(meta_GSE14975$disease)
  # human left atrial appendage - atrial fibrillation patient 
  # 5 
  # human left atrial appendage - control 
  # 5 
  
  res <- purrr::map(search, function(x){
    t.test(expr_GSE14975_final[,x] ~ meta_GSE14975$disease) 
  })
  names(res) <- search
  
  })
}
