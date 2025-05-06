# data type: 基于PEA的蛋白组数据
# company: 汉诺生物
# sample: cell
# software: signature

# Olink 蛋白质组原理介绍
# Olink蛋白质组学的基础是Proximity Extension Assay (PEA)技术，PEA是一种双重识别免疫分析，即两个匹配的抗体同时与一个靶标结合，并以独特的DNA寡核苷酸标记。溶液中的蛋白质使这两种抗体接近，使其 DNA 寡核苷酸杂交，作为 DNA 聚合酶依赖性延伸步骤的模板。这将创建一个双链 DNA"条形码"，该条形码对于特异性抗原来说是唯一的，并且在数量上与靶蛋白的初始浓度成正比。杂交和延伸后立即进行 PCR 扩增 ，最后通过微流控qPCR对扩增子进行定量。
# PEA技术由 Olink Proteomics 公司进行商业化，以开发其生物标志物面板系列。PEA成功地将基于抗体的免疫分析与定量实时PCR（qPCR）和高通量测序（"Next-generation" sequencing technology，NGS）的强大特性相结合，从而形成了一种多重通道和高特异性的方法，目前最多可以同时定量检测3072 种蛋白质生物标志物。PEA是一种非常适合于大规模的精准蛋白质组学研究工具，具有良好的灵敏度、快捷性、高通量分析、以及在多重通道水平上的特异性。

# 项目原始数据
# Olink_**_NPX.xlsx为Olink分析专用软件signature的导出数据，定量数据以蛋白对应荧光信号值经过数据转化并进行log2处理，得到NPX值(Normalization Protein expression，NPX）来呈现。经过log2处理后，部分定量数值可能为负，此为正常情况。经过详细的数据质控可参考质控报告（Summary/Olink_**.pdf)

rm(list = ls());gc()
# load pkgs---------------------------------------------------------------------
require(ggplot2)
require(dplyr)
require(stringr)
# load data---------------------------------------------------------------------
dat <- readxl::read_excel("rawdata.xlsx")

dim(dat)
# [1] 92  6

dat <- dat %>% 
  dplyr::select(-c('Panel','Uniprot ID','OlinkID')) %>% 
  tibble::column_to_rownames("Assay") %>% 
  as.data.frame

group <- dat %>% colnames %>% str_remove("[0-9]") %>% 
  factor(levels = c("BLANK","MODEL","TREAT"))
# QC----------------------------------------------------------------------------
pca <- function(data = dat, group = group){
  library(FactoMineR)##没有请先安装
  pre.pca <- PCA(t(data), graph = FALSE)
  factoextra::fviz_pca_ind(pre.pca,
                           mean.point = F,
                           pointsize = 8,
                           geom= "point",
                           col.ind = group,
                           addEllipses = F) +
    stat_ellipse(level = 0.95, type = "norm", geom = "polygon", col = "black", alpha = 0) +
    # geom_point(size = 4, aes(color = factor(bile_group$Group, labels = c("B","M")))) +
    ggtitle("") +
    # scale_shape_manual(values = c(16,16)) +
    scale_color_manual(values = c("#006eb3","#ca1f25","#1d953f")) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 8, shape = 15)), size = "none", shape = "none") +
    labs(col = "Group") +
    ggtitle("") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          panel.grid = element_blank(), plot.margin = ggplot2::margin(40,40,40,40), ,
          legend.background = element_rect(linetype = 1, color = "black"),
          legend.position = "right",
          axis.text = element_text(size= 16, color = "black"))
}
p1 <- pca(data = dat, group = group)
ggsave(plot = p1, filename = paste0("PCA",".pdf"), 
       width = 7, height = 5.5, dpi = 300)
# DEA---------------------------------------------------------------------------
# BLANK vs MODEL
group <- str_remove(colnames(dat), "[0-9]+$")

report <- data.frame(
  protein = rownames(dat),
  MODEL_mean = rowMeans(dat[,group == "MODEL"]),
  MODEL_sd   = apply(dat[,group == "MODEL"], 1, sd),
  BLANK_mean  = rowMeans(dat[,group == "BLANK"]),
  BLANK_sd    = apply(dat[,group == "BLANK"], 1, sd),
  p.value  = apply(dat, 1, function(x){t.test(x[group == "MODEL"], x[group == "BLANK"])$p.value}),
  stringsAsFactors = FALSE,
  row.names = rownames(dat)
) %>% 
  mutate(FC = MODEL_mean / BLANK_mean,
         fdr = p.adjust(p.value, method = "fdr")) %>% 
  mutate(log2FC = log2(FC))

sum(report$p.value < 0.05) 
sum(report$fdr < 0.05) 
sum(report$fdr < 0.05 & abs(report$log2FC) > log2(1.2)) 



