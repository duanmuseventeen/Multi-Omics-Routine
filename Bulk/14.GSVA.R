# https://blog.csdn.net/dolanjiang/article/details/144673628
# According to the theory, the input should be data adjusted for the length of gene (eg. RPKM, FPKM, TPM, CPM), but not count.
# In fact, the results from TPM or raw count are simliar.

require(GSVA)
require(dplyr)

data <- pseudocounts
dim(data)
# [1] 25644    10

# geneset1 <- qusage::read.gmt("c2.cp.kegg.v7.0.symbols.gmt") #返回的是list
# geneset2 <- clusterProfiler::read.gmt("c2.cp.kegg.v7.0.symbols.gmt") #返回的是表格
# geneset3 <- GSA::GSA.read.gmt("c2.cp.kegg.v7.0.symbols.gmt") #"GSA.genesets", 返回的是list
# geneset4 <- GSEABase::getGmt("c2.cp.kegg.v7.0.symbols.gmt") #list中的list, "GeneSetCollection"
# geneset5 <- cogena::gmt2list("c2.cp.kegg.v7.0.symbols.gmt") #list
# geneset6 <- mogsa::prepMsigDB("c2.cp.kegg.v7.0.symbols.gmt") #list
hallmarks <- qusage::read.gmt("h.all.v2023.1.Hs.symbols.gmt") #返回的是表格

# old version
# gsva_matrix <- GSVA::gsva(as.matrix(count), hallmarks, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

# new version
gsvaPar <- ssgseaParam(exprData = as.matrix(count), 
                       geneSets = hallmarks,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)

# visualization-----------------------------------------------------------------
annotation_cols <- data.frame(
  Group = c(rep("MRP", 4), rep("non-MRP", 6)) %>% 
    factor(levels = c("non-MRP","MRP")),
  row.names = colnames(gsva_data)
)
pheatmap::pheatmap(gsva_data, 
                   clustering_method = "ward.D",
                   annotation_col = annotation_cols,
                   annotation_colors = list(
                     Group = c(MRP = "#C00000", `non-MRP` = "black")),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "row")
