# data type: 基于质谱的蛋白组数据
# company: 景杰生物
# sample: tissue from mice
# software: MaxQuant
# Maxquant处理，LFQ normalization
# 4.MaxQuant中的Intensity，LFQ和iBAQ
# 大佬的软件，三种定量算法都发了文章。
# 
# Intensity是将某Protein Groups里面的所有Unique和Razor peptides的信号强度加起来，作为一个原始强度值。用得很少。
# iBAQ是在Intenstiy的基础上，将原始强度值除以本蛋白的理论肽段数目。一般用于样本内不同蛋白的比较，因为它表征的是蛋白的摩尔比值（copy number）。也可用于不同样本比较，即通过归一化手工校准样本间误差：蛋白IBAQ值除以此样品所有蛋白的强度的和，计算比例（这也是组学中“等质量上样”和“等体积上样”的核心区别，等质量上样来看的是比例，但是计算比例是有压缩效应的）。用得较少。
# LFQ则是将原始强度值在样本之间进行校正，以消除处理、上样、预分、仪器等造成的样本间误差。一般用于同一蛋白不同样本间的比较。不过我们拿到数据后，我们还是会过滤、填充、转换、标准化一条龙走一遍。用得最多。
# 
# 作者：生物信息与育种
# 链接：https://www.jianshu.com/p/de25afe02a33
# 来源：简书
# 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。

rm(list = ls());gc()
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
# Load data---------------------------------------------------------------------
# Normalized by LFQ using MaxQuant
rawdat <- readxl::read_excel("rawdata.xlsx")
# QC----------------------------------------------------------------------------
dat <- rawdat %>% 
  dplyr::select(all_of(c("Gene name",
                         "CKO2_11","CKO2_8","CKO4_19",
                         "WT4_10","WT4_11","WT4_12"))) %>% 
  filter(`Gene name` != "---") %>% 
  filter(complete.cases(CKO2_11)) %>% 
  filter(complete.cases(CKO2_8)) %>% 
  filter(complete.cases(CKO4_19)) %>% 
  filter(complete.cases(WT4_10)) %>% 
  filter(complete.cases(WT4_11)) %>% 
  filter(complete.cases(WT4_12))

# Here, we calculate the mean of the proteins, who have the same name
dat.uni <- dat %>% 
  tidyr::pivot_longer(cols = c("CKO2_11","CKO2_8","CKO4_19","WT4_10","WT4_11","WT4_12"), 
                      names_to = "ID", values_to = "intensity") %>% 
  group_by(`Gene name`, ID) %>% 
  mutate(mean = mean(intensity)) %>% 
  dplyr::select(-intensity) %>% 
  filter(!duplicated(`Gene name`)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = "ID", values_from = "mean") %>% 
  filter(!duplicated(`Gene name`)) %>% 
  filter(complete.cases(`Gene name`)) %>% 
  tibble::column_to_rownames("Gene name")

dim(rawdat)
# [1] 3712   17
dim(dat.uni)
# [1] 3129    6

# Student's t test--------------------------------------------------------------
# In reference to 24942700, we can use Student's t test or Welch test for signifcance.
group <- str_remove(colnames(dat.uni), "[0-9]+_[0-9]+$")

report <- data.frame(
  protein = rownames(dat.uni),
  CKO_mean = rowMeans(dat.uni[,group == "CKO"]),
  CKO_sd   = apply(dat.uni[,group == "CKO"], 1, sd),
  WT_mean  = rowMeans(dat.uni[,group == "WT"]),
  WT_sd    = apply(dat.uni[,group == "WT"], 1, sd),
  p.value  = apply(dat.uni, 1, function(x){t.test(x[group == "CKO"], x[group == "WT"])$p.value}),
  stringsAsFactors = FALSE,
  row.names = rownames(dat.uni)
) %>% 
  mutate(FC = CKO_mean / WT_mean,
         fdr = p.adjust(p.value, method = "fdr")) %>% 
  mutate(log2FC = log2(FC))

sum(report$p.value < 0.05) 
sum(report$fdr < 0.05)
sum(report$fdr < 0.05 & abs(report$log2FC) > log2(1.2))
# Visualization-----------------------------------------------------------------
geom_volcano <- function(dat = report, pos.num = 2,neg.num = -2){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    # tibble::rownames_to_column("protein") %>% 
    mutate(color = case_when(log2FC >  log2(1.2) & fdr < 0.05 ~ 2,
                             log2FC < -log2(1.2) & fdr < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FC))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$protein[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c(NA,0,0,0,0,1,1,NA,NA)
    dat_up$log2FC <- dat_up$log2FC %>% as.numeric
    dat_up$fdr <- dat_up$fdr %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FC)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$protein[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,0,0,1,1,NA,NA)
    dat_down$log2FC <- dat_down$log2FC %>% as.numeric
    dat_down$fdr <- dat_down$fdr %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FC, y = -log10(fdr), col = log2FC, label = protein)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
    )+
    geom_vline(xintercept = c(-log2(1.2),log2(1.2)), color = "gray80", linetype = 2) +
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
      nudge_x = pos.num - as.numeric(dat_up$log2FC),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FC),
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
p <- geom_volcano(dat = report)
ggsave(p, filename = "volcano plot.pdf", 
       width=8, height=8, units="in")

p <- dat.uni[report$protein[report$fdr < 0.05 & abs(report$log2FC) > log2(1.2)],] %>% 
  pheatmap::pheatmap(clustering_method = "ward.D2",
                     scale = "row", 
                     show_rownames = FALSE) %>% 
  ggplotify::as.ggplot()
ggsave(p, filename = "heatmap.pdf", 
       width=6, height=8, units="in")
# Enrichment--------------------------------------------------------------------
# Given the low coverage of protemoics, the GSEA is not suggested.
# ORA---------------------------------------------------------------------------
geneList <- clusterProfiler::bitr(
  report$protein[report$fdr < 0.05 & abs(report$log2FC) > log2(1.2)],
  fromType = "SYMBOL", toType = "ENTREZID", 
  OrgDb = "org.Mm.eg.db")

ORA_KEGG <- enrichKEGG(
  gene = geneList$ENTREZID,
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
  gene = geneList$ENTREZID,
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
dotplot(ORA_GO, showCategory = 20) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

ORA_Reactome <- ReactomePA::enrichPathway(
  gene= geneList$ENTREZID, 
  pvalueCutoff = 1, 
  qvalueCutoff = 1,
  readable=TRUE, 
  organism = "mouse")
ORA_Reactome <- setReadable(ORA_Reactome, 'org.Mm.eg.db', 'ENTREZID')
dotplot(ORA_Reactome, showCategory = 10)

save(ORA_KEGG,ORA_GO,ORA_Reactome, file = "ORA.Rdata")
export::table2excel(ORA_KEGG@result,"ORA_KEGG.xlsx", add.rownames = TRUE)
export::table2excel(ORA_GO@result, "ORA_GO.xlsx", add.rownames = TRUE)
export::table2excel(ORA_Reactome@result, "ORA_REACTOME.xlsx", add.rownames = TRUE)





