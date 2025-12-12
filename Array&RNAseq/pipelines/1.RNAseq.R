setwd("")
rm(list = ls()); gc()
# load pkgs---------------------------------------------------------------------
set.seed(1011)
require(stringr)
require(DESeq2)
require(Mfuzz)
require(dplyr)
require(clusterProfiler)
require(impute)
require(biomaRt)
# load func---------------------------------------------------------------------
geom_volcano <- function(dat,pos.num = 10,neg.num = -10){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    mutate(color = case_when(log2FoldChange >  1 & padj < 0.05 ~ 2,
                             log2FoldChange < -1 & padj < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FoldChange))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$SYMBOL[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c(NA,0,0,0,0,1,1,NA,NA)
    dat_up$log2FoldChange <- dat_up$log2FoldChange %>% as.numeric
    dat_up$padj <- dat_up$padj %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FoldChange)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$SYMBOL[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,0,0,1,1,NA,NA)
    dat_down$log2FoldChange <- dat_down$log2FoldChange %>% as.numeric
    dat_down$padj <- dat_down$padj %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), col = log2FoldChange, label = SYMBOL)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
    )+
    geom_vline(xintercept = c(-1,1), color = "gray80", linetype = 2) +
    geom_hline(yintercept = c(1.30103), color = "gray80", linetype = 2) +
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab("Log2(Fold Change)") +
    labs(color = "") +
    scale_color_gradient2(high = "red3", mid = "white", low = "blue3",
                          midpoint = 0, na.value = "grey80"
                          #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
    )+
    scale_size_continuous(range = c(0.1, 4)) +
    geom_text_repel(
      data = dat_up,
      color = "red3",
      size = 5,
      nudge_x = pos.num - as.numeric(dat_up$log2FoldChange),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FoldChange),
      segment.size = 0.3,
      segment.color = "grey",
      direction="y",
      hjust= 1,
      max.overlaps = Inf) +
    # labs(size = expression("Abundance (log2)"),
    #      color = expression("Direction signed"),
    #      title = trait.names[i]) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
          #legend.spacing.y = unit(1,"cm")
    ) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "italic", colour="grey50", size=20))
}
# load data--------------------------------------------------------------------
rna <- readxl::read_excel("转录组/count.xlsx") %>% 
  as.data.frame
meta <-readxl::read_excel("C:/D/Labmates/齐桂红/转录组分组表.xlsx")

sum(duplicated(rna$id))
sum(duplicated(rna$gene_Dbxref))

dim(pro)
# [1] 9757   16
# 质控 (delete duplicated symbols) ---------------------------------------------
# qc for expression, e.g. filteration of counts < 1000, does not matter significantly, 
# if conducting DEA with DESeq2
sum(duplicated(rna$id))
# [1] 0
sum(duplicated(rna$gene_Dbxref))
# [1] 0
sum(is.na(rna$id))
# [1] 0
sum(rna$id == "")
# [1] 0

count   <- rna %>%
  filter(coding_type == "protein_coding") %>% 
  dplyr::select(-description, - coding_type, -pathway, -pathway_description,
                -GO_ID, -GO_term, -wiki_ID, -wiki_term, -Reactome_ID, -Reactome_term,
                -TF_family) %>% 
  dplyr::select(-gene_Dbxref) %>% 
  tibble::column_to_rownames("id")
# DEA --------------------------------------------------------------------------
group    <- data.frame(
  sample = colnames(count),
  stringsAsFactors = FALSE
) %>% left_join(meta, by = "sample")
group    <- factor(group$group, levels = c("Control","Model","Api"), labels = c("Control","DSS","API")) 
coldata  <- data.frame(row.names = colnames(count), group)
dds      <- DESeqDataSetFromMatrix(countData = count,
                                   colData = coldata,
                                   design = ~group)

dds$group<- relevel(dds$group, ref = "Control") 
dds      <- DESeq(dds)  
DEG_DvsC <- results(dds, name="group_DSS_vs_Control", independentFiltering = FALSE) %>% as.data.frame
DEG_DvsC <- na.omit(DEG_DvsC)

dds$group<- relevel(dds$group, ref = "DSS") 
dds      <- DESeq(dds)  
DEG_AvsD <- results(dds, name="group_API_vs_DSS", independentFiltering = FALSE) %>% as.data.frame
DEG_AvsD <- na.omit(DEG_AvsD)

save(DEG_AvsD, DEG_DvsC, file = "转录组-差异分析.Rdata")
# Visualization ----------------------------------------------------------------
p1 <- geom_volcano(dat = DEG_DvsC, pos.num = 10, neg.num = -10)
p2 <- geom_volcano(dat = DEG_AvsD, pos.num = 10, neg.num = -10)
ggsave(p1, filename = "火山图(DSS vs Control).pdf", width=10, height=6.5, units="in")
ggsave(p2, filename = "火山图(API vs DSS).pdf", width=10, height=6.5, units="in")


DEG <- DEG_DvsC %>% 
  tibble::rownames_to_column("symbol") %>% 
  left_join(DEG_AvsD %>% 
              tibble::rownames_to_column("symbol"),
            by = "symbol", suffix = c(".DvsC", ".AvsD"))

cor.test(DEG$log2FoldChange.AvsD, DEG$log2FoldChange.DvsC)

p3 <- DEG %>% 
  mutate(sig = case_when(
    log2FoldChange.DvsC > 1 & log2FoldChange.AvsD < -1 & padj.DvsC < .05 & padj.AvsD < .05 ~ 2,
    log2FoldChange.DvsC < -1 & log2FoldChange.AvsD > 1 & padj.DvsC < .05 & padj.AvsD < .05 ~ 1,
    TRUE ~ 0 
  )) %>% 
  mutate(sig = factor(sig, labels = c("Non-sig", "Decrease", "Increase"))) %>% 
  ggplot(aes(x = log2FoldChange.DvsC, y = log2FoldChange.AvsD, col = sig)) +
  geom_vline(xintercept = c(-1,1), col = "gray90", linetype = 2) +
  geom_hline(yintercept = c(-1,1), col = "gray90", linetype = 2) +
  geom_point(size = 0.5) +
  annotate(geom = "text", x = -15, y = -5, label = "rho = -0.51\npval < 2.2e-16", hjust = 0) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, col = "steelblue", linewidth = 0.5, linetype = 1) + 
  scale_color_manual(values = c("gray80", "#0C56A0", "#b9292b")) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 16),
        legend.position = "inside", legend.position.inside = c(0.8,0.8))
ggsave(p3, filename = "FC vs FC.pdf", width=8, height=8, units="in")
write.csv(DEG, "DEG.csv")
# Enrichment::ORA + GSEA--------------------------------------------------------
enrichres <- lapply(list(`DSS vs Control` = DEG_DvsC,
                         `API vs DSS` = DEG_AvsD), 
                    function(DEG){
                      DEG2 <- DEG %>% 
                        mutate(id = rownames(DEG)) %>% 
                        left_join(rna %>% 
                                    mutate(entrezgene_id = gene_Dbxref) %>% 
                                    dplyr::select(id, entrezgene_id), by = "id") %>% 
                        dplyr::filter(complete.cases(entrezgene_id)) %>% 
                        mutate(sig = case_when(abs(log2FoldChange) >= 1 & padj < .05 ~ 1,
                                               TRUE ~ 0)) %>% 
                        arrange(desc(log2FoldChange))
                      
                      ORA_KEGG <- enrichKEGG(
                        gene = DEG2$entrezgene_id,
                        organism = "mmu",
                        keyType = "kegg",
                        pvalueCutoff = 1, 
                        qvalueCutoff = 1)
                      ORA_KEGG <- setReadable(ORA_KEGG, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
                      
                      ORA_GO <- enrichGO(
                        gene = DEG2$entrezgene_id,
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
                      
                      ORA_Reactome <- ReactomePA::enrichPathway(
                        gene= DEG2$entrezgene_id, 
                        pvalueCutoff = 1, 
                        readable=TRUE, 
                        organism = "mouse")
                      
                      ORA = list(KEGG = ORA_KEGG, GO = ORA_GO, REACTOME = ORA_Reactome)
                      # GSEA -------------------------------------------------------------------------
                      geneList2 <- DEG2$log2FoldChange
                      names(geneList2) <- DEG2$entrezgene_id
                      geneList2 <- geneList2[!(names(geneList2) %in% names(geneList2)[duplicated(names(geneList2))])]
                      
                      GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                                           seed = 1011,
                                           organism = "mmu",
                                           minGSSize    = 10,
                                           pvalueCutoff = 1,
                                           verbose      = FALSE, 
                                           eps = 0)
                      GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
                      
                      GSEA_GO <- gseGO(geneList = geneList2,
                                       ont = "ALL",
                                       OrgDb = "org.Mm.eg.db",
                                       minGSSize    = 10,
                                       pvalueCutoff = 1,
                                       verbose      = FALSE, 
                                       eps = 0)
                      GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
                      
                      GSEA_REACTOME <- ReactomePA::gsePathway(geneList2,
                                                              organism= "mouse",
                                                              minGSSize= 10,
                                                              maxGSSize= 500,
                                                              pvalueCutoff= 1,
                                                              pAdjustMethod= "BH",
                                                              verbose= FALSE,
                                                              eps= 0)
                      GSEA_REACTOME <- setReadable(GSEA_REACTOME, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
                      
                      GSEA <- list(KEGG = GSEA_KEGG, GO = GSEA_GO, REACTOME = GSEA_REACTOME)
                      
                      return(list(ORA = ORA, GSEA = GSEA))
                    })
# save -------------------------------------------------------------------------
save(enrichres,  file = 'DESeq2-enrichment.Rdata')
write.csv(DEG_DvsC, "DEG_DvsC.csv")
write.csv(DEG_AvsD, "DEG_AvsD.csv")
tag <- c("KEGG", "GO", "REACTOME")
for(i in 1:3){write.csv(enrichres$`DSS vs Control`$ORA [[i]]@result, paste0("DSS vs Control","_ORA",tag[i],".csv"))}
for(i in 1:3){write.csv(enrichres$`DSS vs Control`$GSEA[[i]]@result, paste0("DSS vs Control","_GSEA",tag[i],".csv"))}
for(i in 1:3){write.csv(enrichres$`API vs DSS`$ORA [[i]]@result, paste0("API vs DSS","_ORA",tag[i],".csv"))}
for(i in 1:3){write.csv(enrichres$`API vs DSS`$GSEA[[i]]@result, paste0("API vs DSS","_GSEA",tag[i],".csv"))}
