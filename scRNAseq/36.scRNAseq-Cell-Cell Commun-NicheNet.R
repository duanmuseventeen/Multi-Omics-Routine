# NicheNet
# use data from GSE115469

Learning to use nichenetr
The following vignettes contain the explanation on how to perform a basic NicheNet analysis on a Seurat object. This includes prioritizing ligands and predicting target genes of prioritized ligands. We recommend starting with the step-by-step analysis, but we also demonstrate the use of a single wrapper function. This demo analysis takes only a few minutes to run.

Perform NicheNet analysis starting from a Seurat object: step-by-step analysis:vignette("seurat_steps", package="nichenetr")
Perform NicheNet analysis starting from a Seurat object:vignette("seurat_wrapper", package="nichenetr")
Case study on HNSCC tumor which demonstrates the flexibility of NicheNet. Here, the gene set of interest was determined by the original authors, and the expression data is a matrix rather than a Seurat object.

NicheNet’s ligand activity analysis on a gene set of interest: vignette("ligand_activity_geneset", package="nichenetr")
The following vignettes contain explanation on how to do some follow-up analyses after performing the most basic analysis:
  
  Prioritization of ligands based on expression values: vignette("seurat_steps_prioritization", package="nichenetr")
Inferring ligand-to-target signaling paths: vignette("ligand_target_signaling_path", package="nichenetr")
Assess how well top-ranked ligands can predict a gene set of interest: vignette("target_prediction_evaluation_geneset", package="nichenetr")
Single-cell NicheNet’s ligand activity analysis: vignette("ligand_activity_single_cell", package="nichenetr")
If you want to make a circos plot visualization of the NicheNet output to show active ligand-target links between interacting cells, you can check following vignettes:
  
  Seurat Wrapper + circos visualization:vignette("seurat_wrapper_circos", package="nichenetr").
HNSCC case study + double circos visualization:vignette("circos", package="nichenetr").
People interested in building their own models or benchmarking their own models against NicheNet can read one of the following vignettes:
  
  Model construction: vignette("model_construction", package="nichenetr")
Using LIANA ligand-receptor databases to construct the ligand-target model: vignette("model_construction_with_liana", package="nichenetr")
Model evaluation: target gene and ligand activity prediction: vignette("model_evaluation", package="nichenetr")
Parameter optimization via NSGAII-R: vignette("parameter_optimization", package="nichenetr")
FAQ
Check the FAQ page at FAQ NicheNet: vignette("faq", package="nichenetr")

# Reference---------------------------------------------------------------------
# https://github.com/saeyslab/nichenetr
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper.md

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
# NicheNetR=====================================================================


# Tutorial======================================================================
Perform NicheNet analysis starting from a Seurat object
Robin Browaeys 2023-10-02

In this vignette, you can learn how to perform a basic NicheNet analysis on a 
Seurat (v3-v5) object containing single-cell expression data. Assuming you have 
captured the changes in gene expression resulting from your cell-cell communication 
(CCC) process of interest, a NicheNet analysis can help you to generate hypotheses 
about the CCC process. Specifically, NicheNet can predict 
1) which ligands from the microenvironment or cell population(s) (“sender/niche”)
are most likely to affect target gene expression in an interacting cell population
(“receiver/target”) 
2) which specific target genes are affected by which of these predicted ligands.

The wrapper function we will show consists of the same different steps that are 
discussed in detail in the main NicheNet vignette Perform NicheNet analysis 
starting from a Seurat object: step-by-step analysis.Please make sure you 
understand the different steps described in this vignette before performing a 
real NicheNet analysis on your data. We generally recommend the step-by-step 
analysis as it allows users to adapt specific steps of the pipeline to make them 
more appropriate for their data.

To perform a NicheNet analysis, three features are extracted from the input data:
  the potential ligands, 
  the gene set of interest,
  the background gene set. 
This vignette will extract each feature as described in this flowchart:
  
As example expression data of interacting cells, we will use mouse NICHE-seq 
data to explore intercellular communication in the T cell area in the inguinal 
lymph node before and 72 hours after lymphocytic choriomeningitis virus (LCMV) 
infection (Medaglia et al. 2017). We will focus on CD8 T cells as the receiver 
population, and as this dataset contains two conditions (before and after LCMV 
infection), the differentially expressed genes between these two conditions in 
CD8 T cells will be used as our gene set of interest. We will then prioritize 
which ligands from the microenvironment (sender-agnostic approach) and from 
specific immune cell populations like monocytes, dendritic cells, NK cells, 
B cells, and CD4 T cells (sender-focused approach) can regulate and induce these
observed gene expression changes.

The ligand-target matrix and the Seurat object of the processed NICHE-seq 
single-cell data can be downloaded from Zenodo.

library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)

If you would use and load other packages, we recommend to load these 3 packages
after the others.

Read in NicheNet’s networks
The ligand-target prior model, ligand-receptor network, and weighted integrated networks are needed for this vignette. 
The ligand-target prior model is a matrix describing the potential that a ligand
may regulate a target gene, and it is used to run the ligand activity analysis. 
The ligand-receptor network contains information on potential ligand-receptor 
bindings, and it is used to identify potential ligands. Finally, the weighted 
ligand-receptor network contains weights representing the potential that a 
ligand will bind to a receptor, and it is used for visualization.

organism <- "mouse"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network <- readRDS("lr_network_mouse_21122021.rds")
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final_mouse.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final_mouse.rds")

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
## # A tibble: 6 × 2
##   from          to   
##   <chr>         <chr>
## 1 2300002M23Rik Ddr1 
## 2 2610528A11Rik Gpr15
## 3 9530003J23Rik Itgal
## 4 a             Atrn 
## 5 a             F11r 
## 6 a             Mc1r

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##               2300002M23Rik 2610528A11Rik 9530003J23Rik            a          A2m
## 0610005C13Rik  0.000000e+00  0.000000e+00  1.311297e-05 0.000000e+00 1.390053e-05
## 0610009B22Rik  0.000000e+00  0.000000e+00  1.269301e-05 0.000000e+00 1.345536e-05
## 0610009L18Rik  8.872902e-05  4.977197e-05  2.581909e-04 7.570125e-05 9.802264e-05
## 0610010F05Rik  2.194046e-03  1.111556e-03  3.142374e-03 1.631658e-03 2.585820e-03
## 0610010K14Rik  2.271606e-03  9.360769e-04  3.546140e-03 1.697713e-03 2.632082e-03

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 × 3
##   from          to     weight
##   <chr>         <chr>   <dbl>
## 1 0610010F05Rik App    0.110 
## 2 0610010F05Rik Cat    0.0673
## 3 0610010F05Rik H1f2   0.0660
## 4 0610010F05Rik Lrrc49 0.0829
## 5 0610010F05Rik Nicn1  0.0864
## 6 0610010F05Rik Srpk1  0.123

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 × 3
##   from          to            weight
##   <chr>         <chr>          <dbl>
## 1 0610010K14Rik 0610010K14Rik 0.121 
## 2 0610010K14Rik 2510039O18Rik 0.121 
## 3 0610010K14Rik 2610021A01Rik 0.0256
## 4 0610010K14Rik 9130401M01Rik 0.0263
## 5 0610010K14Rik Alg1          0.127 
## 6 0610010K14Rik Alox12        0.128

Read in the expression data of interacting cells
We processed and aggregated the original dataset by using the Seurat alignment 
pipeline. As we created this object using Seurat v3, it has to be updated with 
UpdateSeuratObject. Note that genes should be named by their official 
mouse/human gene symbol. If your expression data has the older gene symbols, 
you may want to use our alias conversion (convert_alias_to_symbols) 
function to avoid the loss of gene names.

seuratObj <- readRDS("../seuratObj.rds")

# For newer Seurat versions, you may need to run the following
seuratObj <- UpdateSeuratObject(seuratObj)

# Convert gene names
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

if(FALSE){seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513
Visualize which cell populations are present: 
  CD4 T cells (including regulatory T cells), 
  CD8 T cells, 
  B cells, 
  NK cells, 
  dendritic cells (DCs) 
  inflammatory monocytes.

# Note that the number of cells of some cell types is very low and 
# should preferably be higher for a real application
seuratObj@meta.data$celltype %>% table() 
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")

# Visualize the data to see to which condition cells belong. The metadata column 
# that denotes the condition (steady-state or after LCMV infection) is here called 
# ‘aggregate’.

seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141
DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")}

# Perform the NicheNet analysis-------------------------------------------------
In this case study, we want to apply NicheNet to predict which ligands expressed
by the microenvironment (sender-agnostic) and immune cells in the T cell area of 
the lymph node (sender-focused) are most likely to have induced the differential
expression in CD8 T cells after LCMV infection. 
In contrary to NicheNet v1 where we only used the “sender-focused” approach, 
we now recommend users to run both the “sender-agnostic” approach and 
“sender-focused” approach. 
These approaches only affect the list of potential ligands that are considered 
for prioritization. As described in the flowchart above, we do not define any 
sender populations in the ‘sender agnostic’ approach but consider all ligands 
for which its cognate receptor is expressed in the receiver population. The 
sender-focused approach will then filter the list of ligands to ones where the 
ligands are expressed in the sender cell population(s).

As described in the main vignette, the pipeline of a basic NicheNet analysis 
consist of the following steps: 
  * 1. Define a set of potential ligands for both the sender-agnostic and 
    sender-focused approach 
  * 2. Define the gene set of interest: these are the 
    genes in the “receiver/target” cell population that are potentially affected
    by ligands expressed by interacting cells (e.g. genes differentially expressed 
    upon cell-cell interaction) 
  * 3. Define the background genes 
  * 4. Perform NicheNet ligand activity analysis: rank the potential ligands 
    based on the presence of their target genes in the gene set of interest 
    (compared to the background set of genes) 
  * 5. Infer target genes and receptors of top-ranked ligands
All these steps are contained in one of three wrapper functions: 
  nichenet_seuratobj_aggregate, 
  nichenet_seuratobj_cluster_de and 
  nichenet_seuratobj_aggregate_cluster_de. 
These functions differ on how the gene set of interest is calculated, as follows:
  
  Function	Gene set of interest	Background genes
nichenet_seuratobj_aggregate	DE between two conditions of the same cell type	All expressed genes in the cell type
nichenet_seuratobj_cluster_de	DE between two cell types	All expressed genes in both cell types
nichenet_seuratobj_aggregate_cluster_de	DE between two cell types from specific conditions	All expressed genes in both cell types

Note: Cell types should be the identities of the seurat object (check using table(Idents(seuratObj)))

nichenet_seuratobj_aggregate: explain differential expression between two conditions
For the sender-agnostic approach the sender is set to ‘undefined’. The receiver cell population is the ‘CD8 T’ cell population, and the gene set of interest are the genes differentially expressed in CD8 T cells after LCMV infection. Thus, the condition of interest is ‘LCMV’, whereas the reference/steady-state condition is ‘SS’. The column containing condition information is ‘aggregate’. The method to calculate differential expression is the standard Seurat Wilcoxon test. To use other methods, users will have to go through step-by-step analysis. The number of top-ranked ligands that are further used to predict active target genes and construct an active ligand-receptor network is 30 (top_n_ligands). The number of target genes to consider per ligand when performing the target gene inference is 200 (top_n_targets). We only retain ligands and receptors that are expressed in at least a predefined fraction of cells in one cluster (expression_pct, default: 10%).

nichenet_output_agnostic <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  sender = "undefined",
  receiver = "CD8 T", 
  condition_colname = "aggregate",
  condition_oi = "LCMV",
  condition_reference = "SS",
  expression_pct = 0.05,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
# For the sender-focused approach, simply provide one or more sender populations:
  
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
  receiver = "CD8 T", 
  condition_colname = "aggregate",
  condition_oi = "LCMV",
  condition_reference = "SS",
  expression_pct = 0.05,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
Note: It is also possible that you want to consider all cell types present as possible sender cell by defining sender = "all". This also includes the receiver cell type, making that you can look at autocrine signaling as well.

Interpret the NicheNet analysis output
We will investigate the output of the sender-focused approach.

names(nichenet_output)
##  [1] "ligand_activities"                      "top_ligands"                            "top_targets"                           
##  [4] "top_receptors"                          "ligand_target_matrix"                   "ligand_target_heatmap"                 
##  [7] "ligand_target_df"                       "ligand_expression_dotplot"              "ligand_differential_expression_heatmap"
## [10] "ligand_activity_target_heatmap"         "ligand_receptor_matrix"                 "ligand_receptor_heatmap"               
## [13] "ligand_receptor_df"                     "geneset_oi"                             "background_expressed_genes"
Ligand activity analysis results
To see the ranking of ligands based on the predicted ligand activity:
  
nichenet_output$ligand_activities
## # A tibble: 127 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Il27        0.682 0.391         0.316    0.445     1
##  2 Ebi3        0.666 0.264         0.189    0.256     2
##  3 Tnf         0.671 0.205         0.131    0.249     3
##  4 Ptprc       0.660 0.198         0.124    0.168     4
##  5 H2-Eb1      0.656 0.195         0.120    0.182     5
##  6 Vsig10      0.649 0.194         0.119    0.170     6
##  7 H2-M3       0.632 0.192         0.118    0.185     7
##  8 Clcf1       0.637 0.175         0.101    0.162     8
##  9 H2-M2       0.634 0.174         0.0989   0.146    11
## 10 H2-T-ps     0.634 0.174         0.0989   0.146    11
## # ℹ 117 more rows
Ligands are ranked based on the area under the precision-recall curve (AUPR) between a ligand’s target predictions and the observed transcriptional response. Although other metrics like the AUROC and pearson correlation coefficient are also computed, we demonstrated in our validation study that the AUPRwas the most informative measure to define ligand activity (this was the Pearson correlation for v1). The vignette on how we performed the validation can be found at Evaluation of NicheNet’s ligand-target predictions.

To get a list of the top 30 ligands:
  
nichenet_output$top_ligands
 [1] "Il27"    "Ebi3"    "Tnf"     "Ptprc"   "H2-Eb1"  "Vsig10"  "H2-M3"   "Clcf1"   "H2-M2"   "H2-T-ps" "H2-T10"  "H2-T22"  "H2-T24"  "H2-T23"
[15] "H2-K1"   "H2-Q4"   "H2-Q6"   "H2-Q7"   "H2-D1"   "H2-Oa"   "Il18bp"  "Sirpa"   "Cd48"    "App"     "Ccl22"   "Siglech" "Ccl5"    "Siglec1"
[29] "Cd320"   "Adam17"
Below we will show visualizations that are in the output object. In some cases (including this one), not all top ligands that are present in top_ligands will be shown in the plot. The left-out ligands are ligands that don’t have target genes with high enough regulatory potential scores, and therefore did not survive the used cutoffs (in the functions get_weighted_ligand_target_links and prepare_ligand_target_visualization that are run internally). To include them, you can increase the number of target genes considered or be less stringent in the used cutoffs (top_n_targets and cutoff_visualization , respectively). In this case, CCl22 (ranked 25th) is missing from the plots.

# To see which sender cell population expresses which of the top-ranked ligands:
nichenet_output$ligand_expression_dotplot

# As you can see, most op the top-ranked ligands seem to be mainly expressed by dendritic cells and monocytes.
# It could also be interesting to see whether some of these ligands are differentially expressed after LCMV infection.
nichenet_output$ligand_differential_expression_heatmap

# Although this ligand differential expression is not used for prioritization and ranking of the ligands (the ranking is only determined based on enrichment of target genes among DE genes in the receiver, CD8T cells), most of the top-ranked ligands also seem to be upregulated themselves in monocytes after viral infection. This is nice additional “evidence” that these ligands might indeed be important.

# Inferred active ligand-target links
# NicheNet also infers active target genes of these top-ranked ligands, best visualized with the following heatmap showing which top-ranked ligands are predicted to have regulated the expression of which differentially expressed genes:
nichenet_output$ligand_target_heatmap

nichenet_output$ligand_target_heatmap +
  scale_fill_gradient2(low = "whitesmoke",high = "royalblue") +
  xlab("Anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")

# If you want, you can also extract the ligand-target links and their regulatory
# potential scores in matrix or data frame format (e.g. for visualization in 
# other ways or output to a csv file).

nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
##         Adar         B2m        Bst2 Calhm6      Cd274      Cxcl10
## Adam17     0 0.000000000 0.008167279      0 0.06549177 0.011094196
## Cd320      0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Siglec1    0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Ccl5       0 0.000000000 0.000000000      0 0.00000000 0.008424993
## Siglech    0 0.008857572 0.011974948      0 0.01257584 0.008780173
## App        0 0.000000000 0.000000000      0 0.04432138 0.000000000
## Cd48       0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Sirpa      0 0.000000000 0.000000000      0 0.00000000 0.007796006
## Il18bp     0 0.000000000 0.000000000      0 0.00000000 0.007808540
## H2.Oa      0 0.000000000 0.000000000      0 0.00000000 0.008143571
nichenet_output$ligand_target_df # weight column = regulatory potential
## # A tibble: 656 × 3
##    ligand target weight
##    <chr>  <chr>   <dbl>
##  1 Il27   Adar    0.163
##  2 Il27   B2m     0.170
##  3 Il27   Bst2    0.111
##  4 Il27   Calhm6  0.129
##  5 Il27   Cd274   0.111
##  6 Il27   Cxcl10  0.178
##  7 Il27   Cxcr4   0.178
##  8 Il27   Ddx58   0.227
##  9 Il27   Ddx60   0.160
## 10 Il27   Dtx3l   0.150
## # ℹ 646 more rows

# To get a list of the top-predicted target genes of the 30 top-ranked ligands:
nichenet_output$top_targets
##  [1] "Adar"     "B2m"      "Bst2"     "Calhm6"   "Cd274"    "Cxcl10"   "Cxcr4"    "Ddx58"    "Ddx60"    "Dtx3l"    "Eif2ak2"  "Gbp2"     "Gbp3"    
## [14] "Gbp7"     "H2-D1"    "H2-K1"    "H2-M3"    "H2-Q6"    "H2-Q7"    "H2-T10"   "H2-T22"   "H2-T23"   "Ifi203"   "Ifi206"   "Ifi208"   "Ifi209"  
## [27] "Ifi213"   "Ifi35"    "Ifi44"    "Ifih1"    "Ifit1bl1" "Ifit2"    "Ifit3"    "Ifit3b"   "Ifitm3"   "Irf1"     "Irf7"     "Irf9"     "Lgals3bp"
## [40] "Ly6e"     "Mndal"    "Mx1"      "Mx2"      "Nampt"    "Nlrc5"    "Nmi"      "Oas2"     "Oas3"     "Parp12"   "Parp14"   "Parp9"    "Pml"     
## [53] "Psmb8"    "Psmb9"    "Psme1"    "Psme2b"   "Rnf213"   "Samhd1"   "Sp110"    "Stat1"    "Stat2"    "Tap1"     "Tapbp"    "Tnfsf10"  "Trafd1"  
## [66] "Ube2l6"   "Xaf1"     "Ddit4"    "Dhx58"    "Gzmb"     "Isg15"    "Lcp1"     "Oas1a"    "Oas1g"    "Rsad2"    "Zbp1"     "Cd47"     "Ctss"    
## [79] "Trim21"   "Cd69"     "H3f3b"    "Id3"      "Vim"      "Isg20"    "Oasl1"    "Hspa5"    "Ifit1"    "Nt5c3"    "Usp18"    "Basp1"    "Plac8"   
## [92] "Sp100"    "Sp140"    "Ubc"

# You can visualize the expression of these target genes as well (only the top 50 are shown here). Because we only focus on CD8 T cells as receiver cells, we will only show expression in these cells. To emphasize that these target genes are differentially expressed, we split cells up in steady-state cells and cells after response to LCMV infection.
DotPlot(seuratObj %>% subset(idents = "CD8 T"),
        features = nichenet_output$top_targets[1:50] %>%
          rev(), split.by = "aggregate") + coord_flip()

VlnPlot(seuratObj %>% subset(idents = "CD8 T"),
        features = c("Ptprc", "H2-M3", "Cxcl10"), split.by = "aggregate", pt.size = 0, combine = TRUE)

# The display the combined plot of ligand activities, expression, differential 
# expression and target genes of ligands:
nichenet_output$ligand_activity_target_heatmap

# Important: the above figure can be considered as one of the most important summary figures of the NicheNet analysis. Here you can see which ligand-receptor pairs have both high differential expression and ligand activity (=target gene enrichment). These are very interesting predictions as key regulators of your intercellular communication process of interest!
  
# Inferred ligand-receptor interactions for top-ranked ligands
# NicheNet also infers the receiver cell receptors of these top-ranked ligands. You can run following command for a heatmap visualization of the ligand-receptor links:
nichenet_output$ligand_receptor_heatmap

# If you want, you can also extract the ligand-receptor links and their interaction confidence scores in matrix or data frame format (e.g. for visualization in other ways or output to a csv file).

nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
##        H2.T24 H2.T23 H2.T22 H2.T10 H2.T.ps H2.Q7
## Il6ra       0      0      0      0       0     0
## Itgb1       0      0      0      0       0     0
## Notch1      0      0      0      0       0     0
## Ptprc       0      0      0      0       0     0
## Spn         0      0      0      0       0     0
## Cd47        0      0      0      0       0     0
## Cd69        0      0      0      0       0     0
## Ccr7        0      0      0      0       0     0
## Dpp4        0      0      0      0       0     0
## Cd247       0      0      0      0       0     0
nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
## # A tibble: 54 × 3
##    ligand receptor weight
##    <chr>  <chr>     <dbl>
##  1 Adam17 Il6ra     0.447
##  2 Adam17 Itgb1     0.454
##  3 Adam17 Notch1    1.05 
##  4 App    Cd74      0.670
##  5 App    Sorl1     0.922
##  6 Ccl22  Ccr7      0.679
##  7 Ccl22  Dpp4      0.717
##  8 Ccl5   Cxcr3     0.848
##  9 Cd320  Jaml      0.507
## 10 Cd320  Tmem167   0.432
## # ℹ 44 more rows

# To get a list of the receptors of the 30 top-ranked ligands:
nichenet_output$top_receptors
##  [1] "Il6ra"    "Itgb1"    "Notch1"   "Cd74"     "Sorl1"    "Ccr7"     "Dpp4"     "Cxcr3"    "Jaml"     "Tmem167"  "Cd2"      "Il6st"    "Il27ra"  
## [14] "Cd8a"     "Klrd1"    "Cd4"      "Cd247"    "Cd47"     "Ptprc"    "Spn"      "Cd69"     "Tnfrsf1b"

# You can visualize the expression of these as well. Because we only focus on CD8 T cells as receiver cells, we will only show expression in these cells.

DotPlot(seuratObj %>% subset(idents = "CD8 T"),
        features = nichenet_output$top_receptors %>% rev(), split.by = "aggregate") +
  coord_flip()

# If you are interested in checking which geneset (and background set of genes) was used during the ligand activity analysis:
nichenet_output$geneset_oi
##   [1] "Ifi27l2b"      "Irf7"          "Ly6a"          "Stat1"         "Ly6c2"         "Ifit3"         "Ifit1"         "Ly6c1"         "Bst2"         
##  [10] "B2m"           "Rnf213"        "Ifit1bl1"      "Plac8"         "Slfn1"         "Ifi209"        "Isg15"         "Igtp"          "Ifi206"       
##  [19] "Shisa5"        "Ms4a4c"        "H2-K1"         "Zbp1"          "Oasl2"         "Isg20"         "Samhd1"        "Ifi208"        "Ms4a6b"       
##  [28] "Trim30a"       "Usp18"         "Mndal"         "H2-T23"        "Slfn8"         "Gbp2"          "Ifi203"        "Iigp1"         "Tmsb4x"       
##  [37] "H2-T22"        "Rsad2"         "Ly6e"          "Rtp4"          "Ifit3b"        "Zfas1"         "Ifit2"         "Phf11b"        "Xaf1"         
##  [46] "Smchd1"        "Daxx"          "Alb"           "Samd9l"        "Actb"          "Parp9"         "Gbp4"          "Lgals3bp"      "Mx1"          
##  [55] "Ifi213"        "Irgm1"         "2410006H16Rik" "Gbp7"          "Cmpk2"         "Dtx3l"         "Slfn5"         "H2-D1"         "Oasl1"        
##  [64] "Herc6"         "Ifih1"         "Rpsa"          "P2ry13"        "Apoa2"         "Irgm2"         "Tapbp"         "Rps8"          "Stat2"        
##  [73] "Ifi44"         "Phf11c"        "Rpl8"          "Psmb8"         "Gm12250"       "Igfbp4"        "Rplp2-ps1"     "Ddx58"         "Rac2"         
##  [82] "Trafd1"        "Sp100"         "Gbp9"          "Pml"           "Oas2"          "Slfn2"         "Psme1"         "Apoe"          "Gas5"         
##  [91] "H2-Q7"         "Basp1"         "Ms4a4b"        "Rps27a"        "Cd52"          "Znfx1"         "Rpl13"         "Ahsg"          "Oas3"         
## [100] "Nt5c3"         "Rnf114"        "Tap1"          "Rps28"         "Oas1a"         "Rplp0"         "Ddx60"         "Vim"           "Gbp6"         
## [109] "Ifi35"         "Itm2b"         "Ctss"          "Tgtp1"         "Trf"           "Pabpc1"        "H2-Q6"         "Parp14"        "Hspa8"        
## [118] "Tor3a"         "Rpl23"         "Mx2"           "Tmbim6"        "Thy1"          "Ncoa7"         "Dhx58"         "Rps10"         "Rps19"        
## [127] "Psmb9"         "Il2rg"         "Etnk1"         "Irf9"          "Rps3a1"        "Gbp10"         "1600014C10Rik" "Parp12"        "Trim30d"      
## [136] "Eif2ak2"       "Eef1b2"        "Eef2"          "Ncf2"          "Npc2"          "Rps2"          "Rps3"          "Sp110"         "Ube2l6"       
## [145] "Nmi"           "Uba7"          "Psmb10"        "Cxcl10"        "Rpl13a"        "Trim30c"       "Nhp2"          "Tbrg1"         "Jaml"         
## [154] "Usp25"         "Tor1aip2"      "Adar"          "Gzma"          "Gm2000"        "Rps18-ps5"     "Cd53"          "Phf11"         "Hspa5"        
## [163] "Cfl1"          "Crip1"         "Slco3a1"       "Tlr7"          "Trim21"        "Gbp8"          "Rpl10"         "Mycbp2"        "Rps16"        
## [172] "Nlrc5"         "Rplp2"         "Acadl"         "Trim12c"       "Rps4x"         "Irf1"          "Psma2"         "Nme2"          "Tut4"         
## [181] "Apobec3"       "Snord12"       "Phip"          "Gzmb"          "Ifitm3"        "Sp140"         "Dusp2"         "Mrpl30"        "Malat1"       
## [190] "H2-M3"         "Gbp3"          "Tmsb10"        "Dtx1"          "Tmem184b"      "Eef1g"         "Rbl1"          "Epb41l4aos"    "Xpo1"         
## [199] "Rgcc"          "Gm9844"        "Rpl35"         "Rps26"         "Il18bp"        "Sdc3"          "Cxcr4"         "Eif3m"         "Treml2"       
## [208] "Rpl35a"        "Lgals8"        "Pdcd4"         "Arrb2"         "Ubc"           "Clic4"         "H2-T10"        "Rpl10a"        "Lcp1"         
## [217] "Cd274"         "Ddit4"         "Cnn2"          "Nampt"         "Ascc3"         "Ms4a6d"        "Cd47"          "Ogfrl1"        "Snord49b"     
## [226] "Ilrun"         "Calhm6"        "Psme2b"        "Hcst"          "Myh9"          "Rps27"         "Mov10"         "Gm15772"       "Arf4"         
## [235] "Arhgdib"       "Ppib"          "Ubb"           "Trim25"        "Tspo"          "Id3"           "Snord35a"      "Zup1"          "Oas1g"        
## [244] "Ms4a6c"        "Rnf8"          "Casp8"         "Tnfsf10"       "Ptpn7"         "Itk"           "Rps27rt"       "Cd69"          "H3f3b"        
## [253] "Nop10"         "Anxa6"         "Hk1"           "Prkcb"         "Iqgap1"        "Keap1"         "Rpl7"          "Parp10"
nichenet_output$background_expressed_genes %>% length()
## [1] 3476
Results of the sender-agnostic approach
# There is no log-fold change or expression plot because we did not define cell types
nichenet_output_agnostic$ligand_activity_target_heatmap


As you can see in this analysis result, many genes DE in CD8 T cells after LCMV infection are strongly predicted type I interferon targets. The presence of a type I interferon signature in the receiver cell type, but the absence of expression of type I interferons in sender cell types, might indicate that type I interferons are expressed by a different, non-profiled cell type, or at a time point before sampling. The latter could make sense, because there always is a time delay between expression of a ligand-encoding gene and the effect of the ligand on a target/receiver cell (i.e. expression of target genes).

Running multiple NicheNet analyses on different receiver cell populations
In some cases, you might be interested in multiple target/receiver cell populations. You can decide to run this for every cell type separately, or in one line of code as demonstrated here (results are the same). As example, we could have been interested in explaining DE between steady-state and LCMV infection in both CD8 and CD4 T cells.

# To run with  all celltypes in the dataset (only when this would make sense biologically!)
# receiver_celltypes_oi <- seuratObj %>% Idents() %>% unique()

receiver_celltypes_oi <- c("CD4 T", "CD8 T")

nichenet_output <- receiver_celltypes_oi %>% 
  lapply(nichenet_seuratobj_aggregate,
         seurat_obj = seuratObj,
         condition_colname = "aggregate",
         condition_oi = "LCMV",
         condition_reference = "SS",
         sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
         ligand_target_matrix = ligand_target_matrix,
         lr_network = lr_network,
         weighted_networks = weighted_networks)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

names(nichenet_output) <- receiver_celltypes_oi
Check which ligands were top-ranked for both CD8T and CD4T and which ligands were more cell-type specific

common_ligands <- intersect(nichenet_output$`CD4 T`$top_ligands, nichenet_output$`CD8 T`$top_ligands)
print("Common ligands:")
## [1] "Common ligands:"
print(common_ligands)
##  [1] "Ebi3"   "Ptprc"  "H2-M3"  "H2-M2"  "H2-T10" "H2-T22" "H2-T23" "Sirpa"  "H2-K1"  "H2-Q4"  "H2-Q6"  "H2-Q7"  "H2-D1"  "Ccl22"  "Cd48"   "App"   
## [17] "Tgfb1"  "Selplg" "Icam1"  "Btla"   "Cd72"   "B2m"    "Hp"     "Itgb2"

cd4_ligands <- nichenet_output$`CD4 T`$top_ligands %>% setdiff(nichenet_output$`CD8 T`$top_ligands)
cd8_ligands <- nichenet_output$`CD8 T`$top_ligands %>% setdiff(nichenet_output$`CD4 T`$top_ligands)

print("Ligands specifically regulating DE in CD4T:")
## [1] "Ligands specifically regulating DE in CD4T:"
print(cd4_ligands)
## [1] "H2-Eb1"  "H2-Oa"   "Il16"    "Fn1"     "H2-DMb1" "H2-DMb2"

print("Ligands specifically regulating DE in CD8T:")
## [1] "Ligands specifically regulating DE in CD8T:"
print(cd8_ligands)
## [1] "Cxcl10" "Adam17" "Cxcl11" "Tgm2"   "Cxcl9"  "Vcan"
nichenet_seuratobj_cluster_de: explain differential expression between two cell types
Unlike the case above where we applied NicheNet to explain differential expression between two conditions in one cell type, here we try to explain differential expression between two cell populations. DE between cell populations are sometimes (partially) caused by communication with cells in the neighborhood, e.g., the differentiation from a progenitor cell to a differentiated cell might be induced by niche cells. A concrete example is discussed in the paper by Bonnardel et al. (2019): Stellate Cells, Hepatocytes, and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes Colonizing the Liver Macrophage Niche.

However, keep in mind that the comparison that you make should be biologically relevant. as in most cases, differential expression between cell populations will be a result of cell-intrinsic properties (i.e. different cell types have a different gene expression profile) and not of an intercellular communication processes. In such a case, it does not make any sense to use NicheNet.

For demonstration purposes, we will change the Seurat object of the same dataset such that it can be used in this setting.

seuratObj <- SetIdent(seuratObj, value = paste(seuratObj$celltype, seuratObj$aggregate, sep = "_"))
Idents(seuratObj) %>% table()
## .
##   CD8 T_SS   CD4 T_SS    Treg_SS       B_SS      NK_SS    Mono_SS      DC_SS CD8 T_LCMV CD4 T_LCMV     B_LCMV  Treg_LCMV    NK_LCMV  Mono_LCMV 
##        393        601         53         38         37         15          4       1252       1961        344        146         94         75 
##    DC_LCMV 
##         14
Now perform the NicheNet analysis to explain differential expression between the ‘affected’ cell population ‘CD8 T cells after LCMV infection’ and the reference cell population ‘CD8 T cells in steady-state’ by ligands expressed by monocytes and DCs after LCMV infection.

nichenet_output <- nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS",
  receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV", "Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis between two receiver cell clusters"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
Check the top-ranked ligands and their target genes:
  
  nichenet_output$ligand_activity_target_heatmap


Check the expression of the top-ranked ligands:
  
  DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") +
  RotatedAxis()


It could be interesting to check which top-ranked ligands are differentially expressed in monocytes after LCMV infection:
  
  Mono_upregulated_ligands <- FindMarkers(seuratObj, ident.1 = "Mono_LCMV", ident.2 = "Mono_SS") %>% 
  rownames_to_column("gene") %>% filter(avg_log2FC > 0.25 & p_val_adj <= 0.05) %>%
  pull(gene) %>% intersect(nichenet_output$top_ligands)

print("Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-SS and CD8T-LCMV are: ")
## [1] "Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-SS and CD8T-LCMV are: "
print(Mono_upregulated_ligands)
## [1] "B2m"    "H2-D1"  "Cxcl10"
