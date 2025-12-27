setwd("")

rm(list = ls()); gc(); save.image()
# 加载软件----------------------------------------------------------------------
library(stringr)
library(dplyr)
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(enrichplot)
library(GOplot)
require(survminer)
require(survival)
# metascape 在线的
# enrichR 在线的
library(clusterProfiler)
library(pheatmap)
library(data.table)

mytheme <- theme_bw() +
  theme(text = element_text(size = 12), legend.title = element_blank(),
        plot.margin = ggplot2::margin(2,30,30,30), 
        panel.grid = element_blank(), legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4),
        legend.position = c(0.8,0.8))
################################################################################
#--------------------------------  START  --------------------------------------
################################################################################
# 合并数据----------------------------------------------------------------------
samples <- c(
  "C1", "C2",  "C3",  "C4",  "C5", "C6", 
  "M1",  "M2", "M3",  "M4",  "M5",  "M6",
  "E1", "E2",  "E3",  "E4",  "E5", "E6"
)
annot   <- c(
  "peak_name","mz","rt","id","id_zhulab","name","formula","confidence_level",
  "smiles","inchikey","isotope","adduct","total_score","mz_error","rt_error_abs",
  "rt_error_rela","ms2_score","iden_score","iden_type","peak_group_id","base_peak",
  "num_peaks","cons_formula_pred","id_kegg","id_hmdb","id_metacyc",
  "stereo_isomer_id","stereo_isomer_name"
)

uterus_C18_pos <- fread("C18/POS/00_annotation_table/table1_identification.csv")
uterus_C18_neg <- fread("C18/NEG/00_annotation_table/table1_identification.csv")
uterus_HILIC_pos <- fread("HILIC/POS/00_annotation_table/table1_identification.csv")
uterus_HILIC_neg <- fread("HILIC/NEG/00_annotation_table/table1_identification.csv")

all(colnames(uterus_C18_pos) == colnames(uterus_C18_neg))
all(colnames(uterus_C18_pos) == colnames(uterus_HILIC_pos))
all(colnames(uterus_C18_pos) == colnames(uterus_HILIC_neg))

uterus_C18_pos  <- uterus_C18_pos   %>% dplyr::select(all_of(c(annot, samples)))
uterus_C18_neg  <- uterus_C18_neg   %>% dplyr::select(all_of(c(annot, samples)))
uterus_HILIC_pos<- uterus_HILIC_pos %>% dplyr::select(all_of(c(annot, samples)))
uterus_HILIC_neg<- uterus_HILIC_neg %>% dplyr::select(all_of(c(annot, samples)))

uterus_C18_pos  <- uterus_C18_pos   %>% mutate(category = "uterus_C18_pos")
uterus_C18_neg  <- uterus_C18_neg   %>% mutate(category = "uterus_C18_neg")
uterus_HILIC_pos<- uterus_HILIC_pos %>% mutate(category = "uterus_HILIC_pos")
uterus_HILIC_neg<- uterus_HILIC_neg %>% mutate(category = "uterus_HILIC_neg")

lapply(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg), dim)
# 质控(过滤未鉴定离子+过滤低丰度代谢物+过滤缺失+缺失值填补+代谢物去重)----------
# 对于先注释后分析的代谢组, 过滤未鉴定离子--------------------------------------
uterus_C18_pos  <- uterus_C18_pos %>% filter(name != "")
uterus_C18_neg  <- uterus_C18_neg %>% filter(name != "")
uterus_HILIC_pos<- uterus_HILIC_pos %>% filter(name != "")
uterus_HILIC_neg<- uterus_HILIC_neg %>% filter(name != "")

lapply(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg), dim)
# 过滤低丰度代谢物--------------------------------------------------------------
cutoff <- 1000

uterus_C18_pos.expr <- uterus_C18_pos %>% dplyr::select(all_of(samples))
uterus_C18_pos.meta <- uterus_C18_pos %>% dplyr::select(-all_of(samples))
uterus_C18_pos.expr[uterus_C18_pos.expr < cutoff] <- NA

uterus_C18_neg.expr <- uterus_C18_neg %>% dplyr::select(all_of(samples))
uterus_C18_neg.meta <- uterus_C18_neg %>% dplyr::select(-all_of(samples))
uterus_C18_neg.expr[uterus_C18_neg.expr < cutoff] <- NA

uterus_HILIC_pos.expr <- uterus_HILIC_pos %>% dplyr::select(all_of(samples))
uterus_HILIC_pos.meta <- uterus_HILIC_pos %>% dplyr::select(-all_of(samples))
uterus_HILIC_pos.expr[uterus_HILIC_pos.expr < cutoff] <- NA

uterus_HILIC_neg.expr <- uterus_HILIC_neg %>% dplyr::select(all_of(samples))
uterus_HILIC_neg.meta <- uterus_HILIC_neg %>% dplyr::select(-all_of(samples))
uterus_HILIC_neg.expr[uterus_HILIC_neg.expr < cutoff] <- NA

uterus_C18_pos  <- cbind(uterus_C18_pos.meta, uterus_C18_pos.expr)
uterus_C18_neg  <- cbind(uterus_C18_neg.meta, uterus_C18_neg.expr)
uterus_HILIC_pos<- cbind(uterus_HILIC_pos.meta, uterus_HILIC_pos.expr)
uterus_HILIC_neg<- cbind(uterus_HILIC_neg.meta, uterus_HILIC_neg.expr)

lapply(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg), dim)

# 过滤缺失----------------------------------------------------------------------
doqc <- function(dat, sele.col = samples) {
  na.ratio <- apply(dat %>% dplyr::select(all_of(sele.col)), 1, 
                    function(x){sum(is.na(x)) / length(sele.col)})
  return(na.ratio)
}

uterus_C18_pos$na.ratio   <- doqc(uterus_C18_pos)
uterus_C18_neg$na.ratio   <- doqc(uterus_C18_neg)
uterus_HILIC_pos$na.ratio <- doqc(uterus_HILIC_pos)
uterus_HILIC_neg$na.ratio <- doqc(uterus_HILIC_neg)

sum(uterus_C18_pos$na.ratio   < 0.2) / nrow(uterus_C18_pos) # [1] 0.8272727
sum(uterus_C18_neg$na.ratio   < 0.2) / nrow(uterus_C18_neg) # [1] 0.884
sum(uterus_HILIC_pos$na.ratio < 0.2) / nrow(uterus_HILIC_pos) # [1] 0.7447496
sum(uterus_HILIC_neg$na.ratio < 0.2) / nrow(uterus_HILIC_neg) # [1] 0.7468354

uterus_C18_pos   <- uterus_C18_pos   %>% filter(na.ratio < 0.2)
uterus_C18_neg   <- uterus_C18_neg   %>% filter(na.ratio < 0.2)
uterus_HILIC_pos <- uterus_HILIC_pos %>% filter(na.ratio < 0.2)
uterus_HILIC_neg <- uterus_HILIC_neg %>% filter(na.ratio < 0.2)

lapply(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg), dim)
# 填补--------------------------------------------------------------------------
uterus_C18_pos.expr  <- uterus_C18_pos %>% dplyr::select(all_of(samples))
uterus_C18_pos.meta  <- uterus_C18_pos %>% dplyr::select(-all_of(samples))
uterus_C18_pos.expr  <- impute::impute.knn(uterus_C18_pos.expr %>% as.matrix)
 
uterus_C18_neg.expr  <- uterus_C18_neg %>% dplyr::select(all_of(samples))
uterus_C18_neg.meta  <- uterus_C18_neg %>% dplyr::select(-all_of(samples))
uterus_C18_neg.expr  <- impute::impute.knn(uterus_C18_neg.expr %>% as.matrix)

uterus_HILIC_pos.expr<- uterus_HILIC_pos %>% dplyr::select(all_of(samples))
uterus_HILIC_pos.meta<- uterus_HILIC_pos %>% dplyr::select(-all_of(samples))
uterus_HILIC_pos.expr<- impute::impute.knn(uterus_HILIC_pos.expr %>% as.matrix)

uterus_HILIC_neg.expr<- uterus_HILIC_neg %>% dplyr::select(all_of(samples))
uterus_HILIC_neg.meta<- uterus_HILIC_neg %>% dplyr::select(-all_of(samples))
uterus_HILIC_neg.expr<- impute::impute.knn(uterus_HILIC_neg.expr %>% as.matrix)

uterus_C18_pos  <- cbind(uterus_C18_pos.meta, uterus_C18_pos.expr$data)
uterus_C18_neg  <- cbind(uterus_C18_neg.meta, uterus_C18_neg.expr$data)
uterus_HILIC_pos<- cbind(uterus_HILIC_pos.meta, uterus_HILIC_pos.expr$data)
uterus_HILIC_neg<- cbind(uterus_HILIC_neg.meta, uterus_HILIC_neg.expr$data)

lapply(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg), dim)

# 合并数据----------------------------------------------------------------------
uterus_dat <- bind_rows(list(uterus_C18_pos, uterus_C18_neg, uterus_HILIC_pos, uterus_HILIC_neg))
dim(uterus_dat) # [1] 1400   48

sum(duplicated(uterus_dat$name))
# [1] 423
# uterus_dat[uterus_dat$smiles == "C1=C(NC(=O)N=C1)N",]

uterus_dat$median <- apply(uterus_dat %>% dplyr::select(samples), 1, median)
uterus_dat <- uterus_dat %>% 
  arrange(desc(median)) %>% 
  filter(!duplicated(name))

dim(uterus_dat)
# [1] 977  49

# View(uterus_dat)
# 可视化--------------------------------------------------------------------------
uterus_dat <- as.data.frame(uterus_dat)
sum(duplicated(uterus_dat$name))
rownames(uterus_dat) <- uterus_dat$name

pca <- function(dat, group, labels, colors = c("#13679E","#AB3A29","#1E7C4A")){
  library(FactoMineR)##没有请先安装
  pre.pca <- PCA(dat, graph = FALSE)
  factoextra::fviz_pca_ind(pre.pca,
                           mean.point = F, 
                           pointsize = 2,
                           geom= c("point","text"),
                           col.ind = group,
                           addEllipses = T) +
    stat_ellipse(level = 0.95, type = "norm", geom = "polygon", alpha = 0) +
    scale_shape_manual(values = rep(16, length(labels))) +
    scale_fill_manual(values = colors[1:length(labels)]) +
    scale_color_manual(values = colors[1:length(labels)]) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 8, shape = 15)), size = "none", shape = "none") +
    ggtitle("") +
    # coord_fixed(ratio = fix) +
    mytheme
}

p_pca <- pca(dat = uterus_dat[,samples] %>% t, 
             group = str_remove_all(colnames(uterus_dat[,samples]), "[0-9]") %>% 
               factor(levels = c("C","M","E")), 
             labels = c("C","M","E"))

ggsave(p_pca, filename = "uterus_PCA.pdf", width=6, height=6, units="in")
# 差异分析----------------------------------------------------------------------
uterus_dat.scale <- uterus_dat[,samples] / rowSums(uterus_dat[,samples]) * 1000000
uterus_dat.scale <- log2(uterus_dat.scale + 1)
group            <- colnames(uterus_dat.scale) %>% str_remove("\\d+")

uterus_res <- data.frame(
  name   = rownames(uterus_dat.scale),
  mean_C = rowMeans(uterus_dat.scale %>% dplyr::select(starts_with("C"))),
  mean_M = rowMeans(uterus_dat.scale %>% dplyr::select(starts_with("M"))),
  mean_E = rowMeans(uterus_dat.scale %>% dplyr::select(starts_with("E"))),
  sd_C   = apply(uterus_dat.scale %>% dplyr::select(starts_with("C")), 1, sd),
  sd_M   = apply(uterus_dat.scale %>% dplyr::select(starts_with("M")), 1, sd),
  sd_E   = apply(uterus_dat.scale %>% dplyr::select(starts_with("E")), 1, sd),
  pval_M_vs_C = apply(uterus_dat.scale, 1, function(x){t.test(x[group == "C"], x[group == "M"], var.equal = TRUE)$p.value }),
  pval_E_vs_M = apply(uterus_dat.scale, 1, function(x){t.test(x[group == "M"], x[group == "E"], var.equal = TRUE)$p.value }),
  row.names = rownames(uterus_dat.scale)
) %>% mutate(
  log2FC_M_vs_C = mean_M - mean_C,
  log2FC_E_vs_M = mean_E - mean_M,
  fdrq_M_vs_C   = p.adjust(pval_M_vs_C, method = "fdr"),
  fdrq_E_vs_M   = p.adjust(pval_E_vs_M, method = "fdr")
)

# check----
head(uterus_res)

t.test(uterus_dat.scale["Linoleic acid", group == "C"],
       uterus_dat.scale["Linoleic acid", group == "M"])
# stat ----
sum(uterus_res$log2FC_M_vs_C > log2(1.5) & uterus_res$fdrq_M_vs_C < 0.05 &
      uterus_res$log2FC_E_vs_M < -log2(1.5) & uterus_res$fdrq_E_vs_M < 0.05)
sum(uterus_res$log2FC_M_vs_C < -log2(1.5) & uterus_res$fdrq_M_vs_C < 0.05 &
      uterus_res$log2FC_E_vs_M > log2(1.5) & uterus_res$fdrq_E_vs_M < 0.05)

sum(uterus_res$log2FC_M_vs_C > log2(1.5) & uterus_res$pval_M_vs_C < 0.05 & 
      uterus_res$log2FC_E_vs_M < -log2(1.5) & uterus_res$pval_E_vs_M < 0.05) # 73
sum(uterus_res$log2FC_M_vs_C < -log2(1.5) & uterus_res$pval_M_vs_C < 0.05 &
      uterus_res$log2FC_E_vs_M > log2(1.5) & uterus_res$pval_E_vs_M < 0.05) # 2

uterus_res2export <- uterus_res %>% 
  left_join(uterus_dat %>% dplyr::select(id_zhulab, name, formula, smiles, id_kegg, id_hmdb, id_metacyc),
            by = "name")

write.csv(uterus_res2export, "uterus_res.csv")
save(uterus_res2export, file = "uterus_res.Rdata")

# HEATMAP ----
myout.IQR <- function(x){
  y <- x
  x[y < (quantile(y, .25) - 1.5 * IQR(y))] <- (quantile(y, .25) - 1.5 * IQR(y))
  x[y > (quantile(y, .75) + 1.5 * IQR(y))] <- (quantile(y, .75) + 1.5 * IQR(y))
  return(x)
}
my.norm   <- function(x){
  x <- (x - min(x))/(max(x) - min(x))
}

dat.heat <- uterus_dat.scale[uterus_res$name[
  uterus_res$log2FC_M_vs_C > log2(1.5) & uterus_res$pval_M_vs_C < 0.05 &
    uterus_res$log2FC_E_vs_M < -log2(1.5) & uterus_res$pval_E_vs_M < 0.05 |
    uterus_res$log2FC_M_vs_C < -log2(1.5) & uterus_res$pval_M_vs_C < 0.05 &
    uterus_res$log2FC_E_vs_M > log2(1.5) & uterus_res$pval_E_vs_M < 0.05
],]

dat.iqr  <- apply(dat.heat, 1, myout.IQR)
dat.iqr  <- dat.iqr[,hclust(dist(t(dat.iqr)),method = "ward.D2")$order]
p_heat   <- pheatmap::pheatmap(dat.iqr %>% t, 
                             color = colorRampPalette(c("#00bdff","white","#ff2a2a"))(200),
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             clustering_method = "ward.D2",
                             scale = "row") %>% 
  ggplotify::as.ggplot()
ggsave(p_heat, filename = "uterus_heat.pdf", width=20, height=10, units="in")

# VOLCANO PLOT----
p1 <- ggplot(uterus_res2export %>% 
               mutate(sig = case_when(log2FC_M_vs_C > log(1.5) & pval_M_vs_C < 0.05 ~ 2,
                                      log2FC_M_vs_C < -log(1.5) & pval_M_vs_C < 0.05 ~ 0,
                                      TRUE ~ 1) %>% factor), 
             aes(x = log2FC_M_vs_C, y = -log10(pval_M_vs_C), 
                 col = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = "gray80") +
  geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = 2, col = "gray80") +
  geom_point() +
  # ggrepel::geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values = c("#13679E","gray80","#AB3A29")) + 
  ylab(expression("-log10("~italic(P~value)~")")) +
  guides(color = "none") +
  # coord_fixed(ratio = fix) +
  mytheme
ggsave(p1, filename = "uterus_vol_MvC.pdf", width=6, height=6, units="in")

p2 <- ggplot(uterus_res2export %>% 
               mutate(sig = case_when(log2FC_E_vs_M > log(1.5) & pval_E_vs_M < 0.05 ~ 2,
                                      log2FC_E_vs_M < -log(1.5) & pval_E_vs_M < 0.05 ~ 0,
                                      TRUE ~ 1) %>% factor), 
             aes(x = log2FC_E_vs_M, y = -log10(pval_E_vs_M), 
                 col = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = "gray80") +
  geom_vline(xintercept = c(-log(1.5), log(1.5)), linetype = 2, col = "gray80") +
  geom_point() +
  # ggrepel::geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values = c("#13679E","gray80","#AB3A29")) + 
  ylab(expression("-log10("~italic(P~value)~")")) +
  guides(color = "none") +
  # coord_fixed(ratio = fix) +
  mytheme
ggsave(p2, filename = "uterus_vol_EvM.pdf", width=6, height=6, units="in")
