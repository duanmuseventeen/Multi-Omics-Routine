# iTALK
# use data from GSE115469

# Reference---------------------------------------------------------------------
# https://github.com/Coolgenome/iTALK
# https://github.com/Coolgenome/iTALK/blob/master/example/example_code.r
# https://www.jianshu.com/p/f60009e3ba57
# https://blog.csdn.net/qq_42090739/article/details/127317871

# http://www.cellchat.org/
# https://github.com/jinworks/CellChat/
# https://www.jianshu.com/p/30c6e8a24415
# https://www.jianshu.com/p/03f9c2f3ee4f
# https://www.jianshu.com/p/38a9376f5286
# https://github.com/Teichlab/cellphonedb
# https://www.jianshu.com/p/dfe35a2a02d4
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(Matrix)
library(scDblFinder)
## CELL CHAT
# devtools::install_github("sqjin/CellChat")
library(CellChat)
# devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
# devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(iTALK)
# BiocManager::install("SingleCellSignalR")
library(SingleCellSignalR)
# devtools::install_github("arc85/celltalker")
library(celltalker)
# devtools::install_github("mkarikom/RSoptSC")
library(RSoptSC)
# iTALK=========================================================================
dat.italk <- scobj.h.sc@assays$RNA$counts %>% 
  t %>%  as.data.frame %>% 
  mutate(
    cell_type = scobj.h.sc@meta.data$cell_type,
    compare_group = scobj.h.sc@meta.data$orig.ident
    )

# This function loads the count data as a dataframe. It assumes that each line 
# contains gene expression profile of one single cell, and each column contains 
# the one single gene expression profile in different cells. The dataframe should 
# also contain the cell type information with column name 'cell_type'. Group 
# information should also be included as 'compare_group' if users want to call 
# differntial expressed ligand-receptor pairs. Batch information as 'batch' is 
# optional. If included, users may want to use the raw count data for later analysis.
highly_exprs_genes <- rawParse(dat.italk, top_genes = 50, stats = 'mean')

comm_list <- c('growth factor','other','cytokine','checkpoint')

cell_col <- structure(
  c("#F7935A", "#FCB461", "#FED297", "#FFF0DC", "#ABE0E4", "#7FC6D3", "#779FC6"),
  names = unique(dat.italk$cell_type)
  )

par(mfrow=c(1,2))
res<-NULL

for(comm_type in comm_list){
  
  # This function loads the highly expressed genes or differentail expressed genes 
  # as a dataframe. Significant interactions are found through mapping these genes 
  # to our ligand-receptor database.
  # ..data_2	Second dataset used to find ligand-receptor pairs. If set NULL, 
  #   paris will be found within data_1. Otherwise, pairs will be found between 
  #   data_1 and data_2. Default is NULL.
  # ..datatype  Type of data used as input. Options are "mean count" and "DEG"
  # ..comm_type Communication type. Available options are "cytokine", "checkpoint", 
  #   "growth factor", "other"
  # ..database	Database used to find ligand-receptor pairs. If set NULL, 
  #   the build-in database will be used.
  res_cat <- FindLR(highly_exprs_genes, datatype = 'mean count', comm_type = comm_type) %>% 
    arrange(desc(cell_from_mean_exprs * cell_to_mean_exprs)) %>% 
    mutate(arr.col = case_when(
      cell_from == "Endothelia" ~ "#F7935A",
      cell_from == "Cholangiocyte" ~ "#FCB461",
      cell_from == "Mononuclear phagocytes" ~ "#FED297",
      cell_from == "Innate lymphoid cell" ~ "#FFF0DC",
      cell_from == "Hepatocyte" ~ "#ABE0E4",
      cell_from == "Mesenchyme" ~ "#7FC6D3",
      cell_from == "B & Plasma" ~ "#779FC6"))
  
  # res_cat <- res_cat[order(res_cat$cell_from_mean_exprs * res_cat$cell_to_mean_exprs, decreasing=T),]
  
  NetView(res_cat, col = cell_col, vertex.label.cex=1, arrow.width=1, edge.max.width=5)
  
  # link.arr.lwd	
  # line width of the single line link which is put in the center of the belt.
  # link.arr.width	
  # size of the single arrow head link which is put in the center of the belt.
  LRPlot(
    res_cat, 
    datatype = 'mean count', 
    cell_col = cell_col,
    link.arr.col = res_cat$arr.col,
    link.arr.lwd = res_cat$cell_from_mean_exprs,
    link.arr.width = res_cat$cell_to_mean_exprs
    )
  
  title(comm_type)
  
  res <- rbind(res, res_cat)
}
# DEA---------------------------------------------------------------------------
# randomly assign the compare group to each sample
data <- data %>% mutate(compare_group = sample(2,nrow(data),replace=TRUE))

# find DEGenes of regulatory T cells and NK cells between these 2 groups
deg_t <- DEG(data %>% filter(cell_type == 'regulatory_t'), method='Wilcox', contrast=c(2,1))
deg_nk<- DEG(data %>% filter(cell_type == 'cd56_nk'), method='Wilcox', contrast=c(2,1))

# find significant ligand-receptor pairs and do the plotting
par(mfrow=c(1,2))
res<-NULL

function (data_1, data_2 = NULL, datatype, comm_type, database = NULL) 
{
  if (is.null(database)) {
    database <- iTALK:::database
  }
  database <- database[database$Classification == comm_type, 
  ]
  if (datatype == "mean count") {
    gene_list_1 <- data_1
    if (is.null(data_2)) {
      gene_list_2 <- gene_list_1
    }
    else {
      gene_list_2 <- data_2
    }
    
    # here, we can know, when data_2 is used, this function only consider 
    # signal from data1 to data2
    
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in% 
                          gene_list_1$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in% 
                            gene_list_2$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_1 <- database[ind, c("Ligand.ApprovedSymbol", 
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_1[, 
                                                                                           c("gene", "exprs", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_from_mean_exprs = exprs, cell_from = cell_type) %>% 
      left_join(gene_list_2[, c("gene", "exprs", "cell_type")], 
                by = c(Receptor.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_to_mean_exprs = exprs, cell_to = cell_type)
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in% 
                          gene_list_2$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in% 
                            gene_list_1$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_2 <- database[ind, c("Ligand.ApprovedSymbol", 
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_2[, 
                                                                                           c("gene", "exprs", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_from_mean_exprs = exprs, cell_from = cell_type) %>% 
      left_join(gene_list_1[, c("gene", "exprs", "cell_type")], 
                by = c(Receptor.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_to_mean_exprs = exprs, cell_to = cell_type)
    FilterTable <- rbind(FilterTable_1, FilterTable_2)
  }
  else if (datatype == "DEG") {
    gene_list_1 <- data_1
    if (is.null(data_2)) {
      gene_list_2 <- gene_list_1
    }
    else {
      gene_list_2 <- data_2
    }
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in% 
                          gene_list_1$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in% 
                            gene_list_2$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_1 <- database[ind, c("Ligand.ApprovedSymbol", 
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_1[, 
                                                                                           c("gene", "logFC", "q.value", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_from_logFC = logFC, cell_from_q.value = q.value, 
                    cell_from = cell_type) %>% left_join(gene_list_2[, 
                                                                     c("gene", "logFC", "q.value", "cell_type")], by = c(Receptor.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_to_logFC = logFC, cell_to_q.value = q.value, 
                    cell_to = cell_type)
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in% 
                          gene_list_2$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in% 
                            gene_list_1$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_2 <- database[ind, c("Ligand.ApprovedSymbol", 
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_2[, 
                                                                                           c("gene", "logFC", "q.value", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_from_logFC = logFC, cell_from_q.value = q.value, 
                    cell_from = cell_type) %>% left_join(gene_list_1[, 
                                                                     c("gene", "logFC", "q.value", "cell_type")], by = c(Receptor.ApprovedSymbol = "gene")) %>% 
      dplyr::rename(cell_to_logFC = logFC, cell_to_q.value = q.value, 
                    cell_to = cell_type)
    FilterTable <- rbind(FilterTable_1, FilterTable_2)
  }
  else {
    stop("Error: invalid data type")
  }
  FilterTable <- FilterTable[!duplicated(FilterTable), ]
  res <- as.data.frame(FilterTable) %>% dplyr::rename(ligand = Ligand.ApprovedSymbol, 
                                                      receptor = Receptor.ApprovedSymbol)
  if (datatype == "DEG") {
    res <- res[!(res$cell_from_logFC == 1e-04 & res$cell_to_logFC == 
                   1e-04), ]
  }
  res <- res %>% mutate(comm_type = comm_type)
  return(res)
}

for(comm_type in comm_list){
  
  res_cat<-FindLR(deg_t, deg_nk, datatype='DEG',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
  
  #plot by ligand category
  if(nrow(res_cat)==0){
    next
  }else if(nrow(res_cat>=20)){
    LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
  }else{
    LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
  }
  
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  title(comm_type)
  
  res<-rbind(res, res_cat)
}