# Cell Chat
# use data from GSE115469

# Reference---------------------------------------------------------------------
# https://uci-genpals.github.io/signaling/2021/02/23/cellchat.html
# https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

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
# Load data---------------------------------------------------------------------
load("GSE115469-scobj.h.sc.annot.Rdata")
# Cell Chat=====================================================================
#以思维导图方式查看pbm3k.final的结构
# devtools::install_github("pzhaonet/mindr")
library(mindr)
(out <- capture.output(str(scobj.h.sc)))

#### 1.Create cellchat obj------------------------------------------------------
# Create a new CellChat object from a data matrix, Seurat or SingleCellExperiment object
# ..object a normalized (NOT count) data matrix (genes by cells), Seurat or SingleCellExperiment object
cellchat <- createCellChat(object = scobj.h.sc, datatype = "RNA", assay = "RNA", group.by = "cell_type")

cellchat
# An object of class CellChat created from a single dataset 
# 20007 genes.
# 7893 cells. 
# CellChat analysis of single cell RNA-seq data! 

str(cellchat)

levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "cell_type")


groupSize <- cellchat@idents %>% table %>%  as.numeric
groupSize
# [1]  457  113  794 3420 1920   40 1149
#### 2.Load LR DB---------------------------------------------------------------
CellChatDB.human	Ligand-receptor interactions in CellChat database for human
CellChatDB.mouse	Ligand-receptor interactions in CellChat database for mouse
CellChatDB.zebrafish	Ligand-receptor interactions in CellChat database for Zebrafish

# Our database CellChatDB is a manually curated database of literature-supported 
# ligand-receptor interactions in both human and mouse. CellChatDB in mouse contains 
# 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine 
# signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions 
# and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 
# validated molecular interactions, including 61.8% of paracrine/autocrine signaling 
# interactions, 21.7% of extracellular matrix (ECM)-receptor interactions and 16.5% 
# of cell-cell contact interactions.

CellChatDB <- CellChatDB.human

str(CellChatDB)
# interaction、complex、cofactor and geneInfo

showDatabaseCategory(CellChatDB)

# 在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多。

#查看可以选择的侧面，也就是上图左中的三种
unique(CellChatDB$interaction$annotation) 
# [1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 

# select "Secreted Signaling" for subsequent analysis
# Subset CellChatDB databse by only including interactions of interest
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

# set the used database in the object
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use 
#### 3.Preparement--------------------------------------------------------------
# 对表达数据进行预处理，用于细胞间的通信分析。首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
# 如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。

# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat, group.by = "cell_type")

# 相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
# Identify over-expressed signaling genes associated with each cell group
# ..data.use	a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
#* ..idents.use	a subset of cell groups used for analysis
# ..features.name	a char name used for storing the over-expressed signaling genes in 'object@var.features[[features.name]]'
# ..thresh.pc = 0
# ..thresh.fc = 0
# ..thresh.p = 0.05
cellchat <- identifyOverExpressedGenes(cellchat)
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat) 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  

cellchat@LR$LRsig

PPI.human	Human Protein-Protein interactions
PPI.mouse	Mouse Protein-Protein interactions
# Project gene expression data onto a protein-protein interaction network
cellchat <- projectData(cellchat, PPI.human) 
# 找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project
#### 4.Establish Cell Chat Interaction Network----------------------------------
# 通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信。

# 推断配体-受体水平细胞通讯网络（结果储存在@net下面，有一个概率值和对应的pval）
# ⚠️这一步也是CellChat比CellPhoneDB多的一步
# 通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率。

# Compute the communication probability/strength between any interacting cell groups
# 根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
# ..type	methods for computing the average gene expression per cell group. 
#   By default = "triMean", producing fewer but stronger interactions; When setting 
#   ‘type = "truncatedMean"', a value should be assigned to ’trim', producing more 
#   interactions.
# ..trim	the fraction (0 to 0.25) of observations to be trimmed from each end of 
#   x before the mean is computed
# ..population.size	whether consider the proportion of cells in each group across 
#   all sequenced cells. Set population.size = FALSE if analyzing sorting-enriched 
#   single cells, to remove the potential artifact of population size. Set 
#   population.size = TRUE if analyzing unsorted single-cell transcriptomes, with 
#   the reason that abundant cell populations tend to send collectively stronger 
#   signals than the rare cell populations.
# ..LR.use	a subset of ligand-receptor interactions used in inferring communication network
# 
# ..raw.use whether use the raw data (i.e., 'object@data.signaling') or 
#   the projected data (i.e., 'object@data.project'). Set raw.use = FALSE to use 
#   the projected data when analyzing single-cell data with shallow sequencing 
#   depth because the projected data could help to reduce the dropout effects of 
#   signaling genes, in particular for possible zero expression of subunits of 
#   ligands/receptors.
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# If you really worry about that, I think filter the cluster with few cells before creation is better
function (object, min.cells = 10) 
{
  net <- object@net
  cell.excludes <- which(as.numeric(table(object@idents)) <= 
                           min.cells)
  if (length(cell.excludes) > 0) {
    cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", 
        levels(object@idents)[cell.excludes], "\n")
    net$prob[cell.excludes, , ] <- 0
    net$prob[, cell.excludes, ] <- 0
    object@net <- net
  }
  return(object)
}
# cellchat <- filterCommunication(cellchat, min.cells = 10)

subsetCommunication(cellchat)

# 推断信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的pval）
# 我们可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络。
# Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
cellchat <- computeCommunProbPathway(cellchat)

subsetCommunication(cellchat, slot.name = "netP")
#### 5.Visualization------------------------------------------------------------
# Calculate the aggregated network by counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat)

# Return an updated CellChat object:
# 
# 'object@net$count' is a matrix: rows and columns are sources and targets respectively, and elements are the number of interactions between any two cell groups. USER can convert a matrix to a data frame using the function 'reshape2::melt()'
# 'object@net$weight' is also a matrix containing the interaction weights between any two cell groups
# 'object@net$sum' is deprecated. Use 'object@net$weight'
cellchat@net$count
  B & Plasma Cholangiocyte Endothelia Hepatocyte
B & Plasma                     15            35         60         17
Cholangiocyte                  56           122        145        120
Endothelia                     81           162        169        144
Hepatocyte                     36           131        150        121
Innate lymphoid cell           38            87        115         55
Mesenchyme                    104           171        171        164
Mononuclear phagocytes         35           110        132         78
Innate lymphoid cell Mesenchyme Mononuclear phagocytes
B & Plasma                               31         37                     38
Cholangiocyte                            93        119                    131
Endothelia                              133        153                    163
Hepatocyte                               54        115                     98
Innate lymphoid cell                     61         77                     92
Mesenchyme                              111        148                    147
Mononuclear phagocytes                   82         97                    108
  
cellchat@net$weight
B & Plasma Cholangiocyte  Endothelia  Hepatocyte
B & Plasma             0.0011376787  0.0001307159 0.002146987 0.001022693
Cholangiocyte          0.0003342979  0.0001712926 0.001602424 0.003677983
Endothelia             0.0029896072  0.0009674425 0.016301269 0.013702070
Hepatocyte             0.0060262472  0.0079390517 0.044463353 0.208204407
Innate lymphoid cell   0.0056613101  0.0009865568 0.023767729 0.014739574
Mesenchyme             0.0002750721  0.0003964093 0.001399818 0.008418919
Mononuclear phagocytes 0.0027942403  0.0006941276 0.019083999 0.008489463
Innate lymphoid cell   Mesenchyme Mononuclear phagocytes
B & Plasma                      0.011651297 4.140999e-05            0.005179677
Cholangiocyte                   0.002969641 9.299293e-05            0.002266532
Endothelia                      0.040054798 4.466691e-04            0.015323173
Hepatocyte                      0.045924844 2.278486e-03            0.055863518
Innate lymphoid cell            0.128741203 6.612359e-04            0.028541207
Mesenchyme                      0.001849818 9.164666e-05            0.001780003
Mononuclear phagocytes          0.041784154 3.359943e-04            0.029400111

#计算每种细胞各有多少个
groupSize <- cellchat@idents %>% table %>% as.numeric

# Visualize the inferred cell-cell communication network
# ..vertex.weight	The weight of vertex: either a scale value or a vector
#   Default is a scale value being 1, indicating all vertex is plotted in the same 
#   size;
#   Set 'vertex.weight' as a vector to plot vertex in different size; setting 
#   'vertex.weight = NULL' will have vertex with different size that are portional 
#   to the number of cells in each cell group.

par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(
  cellchat@net$count, 
  vertex.weight = 1, 
  weight.scale = T,
  label.edge= F, # label the number
  title.name = "Number of interactions")
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = 1, 
  weight.scale = T,
  label.edge= F, 
  title.name = "Interaction weights/strength")

netVisual_circle(
  cellchat@net$count, 
  vertex.weight = groupSize, 
  weight.scale = T,
  label.edge= F, 
  title.name = "Number of interactions")
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = groupSize, 
  weight.scale = T,
  label.edge= F, 
  title.name = "Interaction weights/strength")

# count
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    arrow.width = 0.2,
    arrow.size = 0.1, 
    edge.weight.max = max(mat), 
    title.name = rownames(mat)[i])
}

# weight
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    arrow.width = 0.2,
    arrow.size = 0.1, 
    edge.weight.max = max(mat), 
    title.name = rownames(mat)[i])
}

cellchat@netP$pathways

# define a numeric vector （淋系细胞）giving the index of the celltype as targets
vertex.receiver = c(1,2,4,6) 

# B & Plasma          Cholangiocyte             Endothelia 
# 457                    113                    794 
# Hepatocyte   Innate lymphoid cell             Mesenchyme 
# 3420                   1920                     40 
# Mononuclear phagocytes 
# 1149 

# Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
netVisual_aggregate(
  cellchat, 
  signaling = "MHC-I",  # a signaling pathway name
  # color.use	= , the character vector defining the color of each cell group
  # top	= , the fraction of interactions to show
  layout = c("hierarchy", "circle", "chord", "spatial")[1], 
  vertex.receiver = c(3, 7) # when layout hierarchy, vertex.receiver must be a vector with more than 1 element
)

netVisual_aggregate(
  cellchat, 
  signaling = "MHC-I",  # a signaling pathway name
  # color.use	= , the character vector defining the color of each cell group
  # top	= , the fraction of interactions to show
  layout = c("hierarchy", "circle", "chord", "spatial")[2],
  vertex.receiver = NULL
    )

netVisual_aggregate(
  cellchat, 
  signaling = "MHC-I",  # a signaling pathway name
  # color.use	= , the character vector defining the color of each cell group
  # top	= , the fraction of interactions to show
  layout = c("hierarchy", "circle", "chord", "spatial")[3],
  vertex.receiver = NULL
)

netVisual_heatmap(cellchat, signaling = "MHC-I")

# 配体-受体层级的可视化（计算各个ligand-receptor pair对信号通路的贡献）
# 计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对MHC-I贡献最大）
# Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways

netAnalysis_contribution(cellchat) # error
netAnalysis_contribution(cellchat, signaling = "MHC-I") # ok
netAnalysis_contribution(cellchat, signaling = c("MHC-I","FN1")) # ok

# Identify all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway
# ..geneLR.return	whether return the related signaling genes of enriched L-R pairs
extractEnrichedLR(
  cellchat, 
  signaling = "MHC-I", 
  geneLR.return = TRUE,
  enriched.only = TRUE,
  thresh = 0.05)

extractEnrichedLR(
  cellchat, 
  signaling = c("MHC-I","FN1"), 
  geneLR.return = TRUE,
  enriched.only = TRUE,
  thresh = 0.05)

# Registered S3 method overwritten by 'spatstat':#>   method     from#>   print.boxx cli#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.
plotGeneExpression(cellchat, signaling = "MHC-I")
plotGeneExpression(cellchat, signaling = "MHC-I", enriched.only = FALSE)

netVisual_bubble(
  cellchat, 
  sources.use = c(4),
  targets.use = c(1:7), 
  # signaling = c("CCL","CXCL"),
  remove.isolate = FALSE)
#### 6.Systems analysis of cell-cell communication network----------------------
# Compute and visualize the network centrality scores
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-I"))

gg1 + gg2

# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I"))

#### 6.Save the CellChat object-------------------------------------------------
saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")




















