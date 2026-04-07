rm(list = ls())

library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

drug.screen <- readxl::read_excel("drug screen score.xlsx")
reference   <- readxl::read_excel("positive drugs.xlsx")

drug.screen <- drug.screen %>% arrange(desc(score))
drug_list <- drug.screen$score
names(drug_list) <- drug.screen$drug

custom_pathways <- list(Positive_Drugs = reference$`positive drug`)
# fgsea ------------------------------------------------------------------------
gsea_results <- fgsea(pathways = custom_pathways, 
                      stats = drug_list,
                      minSize = 1,      # 即使只有几个阳性药物也能计算
                      maxSize = 500)    

# 查看结果 (NES > 0 且 P < 0.05 表示阳性药物显著富集在打分靠前位置)
print(gsea_results)

plotEnrichment(custom_pathways[["Positive_Drugs"]], drug_list) +
  labs(title = "Enrichment of Positive Drugs",
       subtitle = "Calculated via fgsea") +
  theme_bw()
# GSEA -------------------------------------------------------------------------
custom_geneset <- reference %>% 
  transmute(term = "postive drug", gene = `positive drug`) 

gsea_results_formal <- GSEA(drug_list, 
                            TERM2GENE = custom_geneset, 
                            pvalueCutoff = 1)

gseaplot2(gsea_results_formal, 
          geneSetID = 1, 
          title = "Enrichment of Resistance Module",
          pvalue_table = TRUE)
