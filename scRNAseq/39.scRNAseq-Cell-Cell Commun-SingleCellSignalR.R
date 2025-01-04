# SingleCellSignalR
# use data from GSE115469

SingleCellSignalR
User’s Guide
Simon Cabello-Aguilar1, Jacques Colinge1
1 Institut de Recherche en Cancérologie de Montpellier, Inserm, Montpellier, France ; Institut régional du Cancer Montpellier, Montpellier, France ; Université de Montpellier, Montpellier, France

Introduction
This guide provides an overview of the SingleCellSignalR package, 
a comprehensive framework to obtain cellular network maps from scRNA-seq data. 
SingleCellSignalR comes with a complete pipeline integrating existing methods to 
cluster individual cell transcriptomes and identify cell subpopulations as well as 
novel cellular network-specific algorithms. More advanced users can substitute 
their own logic or alternative tools at various stages of data processing. 
SingleCellSignalR also maps cell subpopulation internal network linked to genes 
of interest through the integration of regulated KEGG and Reactome pathways 
together with ligands and receptors involved in inferred cell-cell interactions. 
The cellular networks can be exported in text files and graphML objects to be further 
explored with Cytoscape (www.cytoscape.org), yEd (www.yworks.com), or 
similar software tools.

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
# ==============================================================================
# ...SingleCellSignalR can be used to generate cellular interaction lists using the 
# cell_signaling() function:

data <- scobj.h.sc@assays$RNA$data

# Computes "autocrine" or "paracrine" interactions between cell clusters.
signal <- cell_signaling(
  # a data frame of n rows (genes) and m columns (cells) of read or UMI counts (note : rownames(data)=genes)
  # log-normalized data, not counts
  data = data, 	
  genes = rownames(data), # a character vector of HUGO official gene symbols of length n
  cluster = scobj.h.sc@meta.data$cell_type %>% factor %>% as.numeric, # a numeric vector of length m
  write = FALSE)

# No such file as table_dge_cluster 1.txt in the cluster-analysis folder
# No such file as table_dge_cluster 2.txt in the cluster-analysis folder
# No such file as table_dge_cluster 3.txt in the cluster-analysis folder
# No such file as table_dge_cluster 4.txt in the cluster-analysis folder
# No such file as table_dge_cluster 5.txt in the cluster-analysis folder
# No such file as table_dge_cluster 6.txt in the cluster-analysis folder
# No such file as table_dge_cluster 7.txt in the cluster-analysis folder
# Paracrine signaling: 
#   Checking for signaling between cell types
# 48 interactions from cluster 1 to cluster 2
# 31 interactions from cluster 1 to cluster 3
# 14 interactions from cluster 1 to cluster 4
# 26 interactions from cluster 1 to cluster 5
# 49 interactions from cluster 1 to cluster 6
# 21 interactions from cluster 1 to cluster 7
# 39 interactions from cluster 2 to cluster 1
# 56 interactions from cluster 2 to cluster 3
# 16 interactions from cluster 2 to cluster 4
# 79 interactions from cluster 2 to cluster 5
# 58 interactions from cluster 2 to cluster 6
# 66 interactions from cluster 2 to cluster 7
# 7 interactions from cluster 3 to cluster 1
# 28 interactions from cluster 3 to cluster 2
# 7 interactions from cluster 3 to cluster 4
# 19 interactions from cluster 3 to cluster 5
# 16 interactions from cluster 3 to cluster 6
# 17 interactions from cluster 3 to cluster 7
# 9 interactions from cluster 4 to cluster 1
# 18 interactions from cluster 4 to cluster 2
# 16 interactions from cluster 4 to cluster 3
# 12 interactions from cluster 4 to cluster 5
# 9 interactions from cluster 4 to cluster 6
# 10 interactions from cluster 4 to cluster 7
# 3 interactions from cluster 5 to cluster 1
# 22 interactions from cluster 5 to cluster 2
# 6 interactions from cluster 5 to cluster 3
# 1 interactions from cluster 5 to cluster 4
# 14 interactions from cluster 5 to cluster 6
# 7 interactions from cluster 5 to cluster 7
# 80 interactions from cluster 6 to cluster 1
# 106 interactions from cluster 6 to cluster 2
# 95 interactions from cluster 6 to cluster 3
# 56 interactions from cluster 6 to cluster 4
# 99 interactions from cluster 6 to cluster 5
# 112 interactions from cluster 6 to cluster 7
# 7 interactions from cluster 7 to cluster 1
# 39 interactions from cluster 7 to cluster 2
# 7 interactions from cluster 7 to cluster 3
# 3 interactions from cluster 7 to cluster 4
# 13 interactions from cluster 7 to cluster 5
# 25 interactions from cluster 7 to cluster 6

class(signal)
[1] "list"
length(signal)
[1] 42
head(signal$`cluster 1-cluster 2`)
cluster 1 cluster 2 interaction type   LRscore
2628      SAA1    ADRA2A        paracrine 0.8714972
2591    RPS27A     ERBB2        paracrine 0.8695930
2959     UBA52     ERBB2        paracrine 0.8632075
1582  HSP90AA1      CFTR        paracrine 0.8613111
1588   HSP90B1     ERBB2        paracrine 0.8598953
1528     HLA-A     ERBB2        paracrine 0.8480552

# An intercellular network can also be generated to map the overall ligand/receptor 
# interactions invoking the inter_network() function:
inter.net <- inter_network(
  data = data, 
  signal = signal, 
  genes = rownames(data), 
  cluster = scobj.h.sc@meta.data$cell_type %>% factor %>% as.numeric, 
  write = FALSE)
No such file as table_dge_cluster 1.txt in the cluster-analysis folder
No such file as table_dge_cluster 2.txt in the cluster-analysis folder
No such file as table_dge_cluster 3.txt in the cluster-analysis folder
No such file as table_dge_cluster 4.txt in the cluster-analysis folder
No such file as table_dge_cluster 5.txt in the cluster-analysis folder
No such file as table_dge_cluster 6.txt in the cluster-analysis folder
No such file as table_dge_cluster 7.txt in the cluster-analysis folder
Paracrine signaling: 
  Checking for signaling between cell types
48 interactions from cluster 1 to cluster 2
31 interactions from cluster 1 to cluster 3
14 interactions from cluster 1 to cluster 4
26 interactions from cluster 1 to cluster 5
49 interactions from cluster 1 to cluster 6
21 interactions from cluster 1 to cluster 7
39 interactions from cluster 2 to cluster 1
56 interactions from cluster 2 to cluster 3
16 interactions from cluster 2 to cluster 4
79 interactions from cluster 2 to cluster 5
58 interactions from cluster 2 to cluster 6
66 interactions from cluster 2 to cluster 7
7 interactions from cluster 3 to cluster 1
28 interactions from cluster 3 to cluster 2
7 interactions from cluster 3 to cluster 4
19 interactions from cluster 3 to cluster 5
16 interactions from cluster 3 to cluster 6
17 interactions from cluster 3 to cluster 7
9 interactions from cluster 4 to cluster 1
18 interactions from cluster 4 to cluster 2
16 interactions from cluster 4 to cluster 3
12 interactions from cluster 4 to cluster 5
9 interactions from cluster 4 to cluster 6
10 interactions from cluster 4 to cluster 7
3 interactions from cluster 5 to cluster 1
22 interactions from cluster 5 to cluster 2
6 interactions from cluster 5 to cluster 3
1 interactions from cluster 5 to cluster 4
14 interactions from cluster 5 to cluster 6
7 interactions from cluster 5 to cluster 7
80 interactions from cluster 6 to cluster 1
106 interactions from cluster 6 to cluster 2
95 interactions from cluster 6 to cluster 3
56 interactions from cluster 6 to cluster 4
99 interactions from cluster 6 to cluster 5
112 interactions from cluster 6 to cluster 7
7 interactions from cluster 7 to cluster 1
39 interactions from cluster 7 to cluster 2
7 interactions from cluster 7 to cluster 3
3 interactions from cluster 7 to cluster 4
13 interactions from cluster 7 to cluster 5
25 interactions from cluster 7 to cluster 6

# At this point the intercellular network have been generated and exported in text 
# and graphML formats in the networks folder.
# A summary of the interactions between cell clusters can be output in the form of 
# a chord diagram by the visualize_interactions() function:
visualize_interactions(signal = signal)

# The details of the interactions between two clusters, for example cluster 1 and 2, 
# can also be shown in the plot window with the visualize_interactions() function. 
# Note that in the example below we ask for the display of two pairs of cell clusters, 
# pair 1 that contains interactions from cluster 1 to 2, and pair 4 from cluster 2 to 1. 
# (names(signal) returns the cell cluster names in each pair, see function 
# visualize_interactions() details.)
visualize_interactions(signal = signal, show.in=c(1,4))

# And these plots can be saved into pdf files in the images folder using the write.in 
# argument of the visualize_interactions() function.
visualize_interactions(signal = signal, write.in=c(1,4))

# ==============================================================================
# SingleCellSignalR package functions have many arguments parameters that 
# can be changed by the user to fit her needs (see Reference Manual for more details). 
# Furthermore, several handy functions that were not illustrated above are 
# provided to generate additional plots or reports.

# Exploiting the cell_classifier clustering
# After running the example in the Quick Start section, the user can define 
# cell clusters after the output of the cell_classifier(). The demo data set
# is comprised of a subset of the 10x PBMC dataset [3], i.e. immune cells. 
# The t-SNE map calculated with the clustering() function will also be used. 
# For this example we will set the plot.details argument to TRUE to monitor 
# the choice of the threshold of gene signature scores.

# Let us use the cell clustering obtained with the cell_classifier() function. 
# Although “undefined” cells may be interesting in some cases, here they form a 
# heterogeneous cluster because they represent cells that seem to be in a 
# transition between two states (“T-cells” and “Cytotoxic cells”, or 
# “Neutrophils” and “Macrophages”, see heatmap above). We discard these cells.

data <- scobj.h.sc@assays$RNA$data
umap <- scobj.h.sc@reductions$UMAP@cell.embeddings
cluster <- scobj.h.sc@meta.data$cell_type %>% factor %>% as.numeric
names(cluster) <- scobj.h.sc@meta.data$cell_type
c.names <- names(cluster) %>% unique

# clust.ana <- cluster_analysis(
#   data = data, 
#   genes = rownames(data), 
#   cluster = cluster, 
#   c.names = c.names,
#   write = FALSE)

# Once the cluster analysis is done, the cell_signaling(), inter_network() 
# functions can be used.
signal <- cell_signaling(
  data = data, 
  genes = rownames(data), 
  cluster = cluster, 
  c.names = c.names,
  write = FALSE)

# We can be interested in genes participating in pathways with a receptor of 
# interest inside a cluster of interest. Let us say ASGR1 in “Macrophages”.
intra <- my_intra_network(
  goi = "PTGDR",	# gene of interest (typically a receptor)
  data = data,
  genes = rownames(data),
  cluster = cluster, 
  coi = "Hepatocyte", 
  c.names = c.names, 
  signal = signal,
  write=FALSE)

# Patwhay(s) that include PTGDR:
#   No associated genes downstream PTGDR in Hepatocyte

# Now, let us take an overview of the signaling between the cell types.
visualize_interactions(signal)

visualize_interactions(signal, show.in = c(3,4))

visualize_interactions(signal, write.in= c(3,4))

# Revise my_intra_network at line ~386 in function RowMeans
# add arg drop = FALSE
my_intra_network <-
function (goi, data, genes, cluster, coi, cell.prop = 0.2, c.names = NULL, 
          signal = NULL, write = TRUE, plot = TRUE, add.lig = TRUE, 
          species = c("homo sapiens", "mus musculus"), connected = FALSE, 
          verbose = TRUE) {
  require(foreach)
  if (dir.exists("networks") == FALSE & write == TRUE) {
    dir.create("networks")
  }
  if (is.null(c.names) == TRUE) {
    c.names <- paste("cluster", seq_len(max(cluster)))
  }
  if (min(cluster) != 1) {
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names) != max(cluster) | sum(duplicated(c.names)) > 
      0 | grepl("/", paste(c.names, collapse = ""))) {
    stop("The length of c.names must be equal to the number of clusters and must\n        contain no duplicates. The cluster names must not include special\n        characters")
  }
  if (!is.element(coi, c.names)) {
    stop(paste(coi, "must be included in c.names.", "If c.names is not provided, it is set to cluster 1, cluster 2, ...,\n        cluster N. WIth N the maximum number of clusters"))
  }
  opar <- par()
  species <- match.arg(species)
  if (species == "mus musculus") {
    Hs2mm <- mm2Hs[, 1]
    mm2Hs <- mm2Hs[, 2]
    names(mm2Hs) <- Hs2mm
    names(Hs2mm) <- as.character(mm2Hs)
    m.names <- mm2Hs[rownames(data)]
    data <- subset(data, (!is.na(m.names)))
    m.names <- m.names[!is.na(m.names)]
    rownames(data) <- as.character(m.names)
    goi <- mm2Hs[goi]
  }
  pw.names <- strsplit(PwC_ReactomeKEGG$pathway, ";")
  pw.sizes <- table(unlist(pw.names))
  max.pw.size <- 500
  good.pw <- pw.sizes[pw.sizes <= max.pw.size]
  data.tmp <- data[, cluster == which(c.names == coi)]
  data.tmp <- data.tmp[rowSums(data.tmp) > 0, ]
  good <- apply(data.tmp, 1, function(x) sum(x > 0)/ncol(data.tmp) > 
                  cell.prop)
  visible.genes <- unique(c(rownames(data.tmp)[good], goi))
  visible.n <- PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn %in% 
                                  visible.genes & PwC_ReactomeKEGG$b.gn %in% visible.genes, 
  ]
  red.visible.n <- simplify_interactions(visible.n, LRdb)
  res <- list()
  qq <- 0
  for (receptors in goi) {
    qq <- qq + 1
    if (!is.element(receptors, rownames(data.tmp))) {
      cat(paste0(receptors, " is not expressed in ", coi), 
          fill = TRUE)
    }
    else {
      contains.receptors <- intersect(unlist(pw.names[PwC_ReactomeKEGG$a.gn %in% 
                                                        receptors | PwC_ReactomeKEGG$b.gn %in% receptors]), 
                                      names(good.pw))
      if (verbose == TRUE) {
        cat(paste0("Patwhay(s) that include ", receptors, 
                   ":"), fill = TRUE)
        for (i in contains.receptors) {
          cat(paste0("   - ", i), fill = TRUE)
        }
      }
      if (sum(grepl("added", contains.receptors)) > 1) {
        contains.receptors <- contains.receptors[contains.receptors != 
                                                   "added"]
      }
      contain.n = NULL
      for (i in contains.receptors) {
        contain.n <- rbind(contain.n, PwC_ReactomeKEGG[grepl(i, 
                                                             PwC_ReactomeKEGG$pathway), ])
      }
      contain.n <- unique(rbind(contain.n, PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn %in% 
                                                              visible.genes & PwC_ReactomeKEGG$b.gn %in% receptors, 
      ], PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn %in% 
                            receptors & PwC_ReactomeKEGG$b.gn %in% visible.genes, 
      ]))
      red.contain.n <- simplify_interactions(contain.n)
      key.visible <- paste(red.visible.n$a.gn, red.visible.n$b.gn, 
                           sep = "|")
      key.contain <- paste(red.contain.n$a.gn, red.contain.n$b.gn, 
                           sep = "|")
      net.n <- red.visible.n[key.visible %in% key.contain, 
      ]
      if (nrow(net.n) > 0) {
        add.net <- NULL
        nam = NULL
        siz = NULL
        if (is.null(signal) == FALSE & add.lig == TRUE) {
          for (i in names(signal)) {
            if (grepl(paste0(coi, "$"), i)) {
              tmp = signal[[i]]
              if (species == "mus musculus") {
                m.1 <- mm2Hs[tmp[, 1]]
                m.2 <- mm2Hs[tmp[, 2]]
                tmp <- subset(tmp, (!is.na(m.1)))
                tmp <- subset(tmp, (!is.na(m.2)))
                m.1 <- m.1[!is.na(m.1)]
                m.2 <- m.2[!is.na(m.2)]
                tmp[, 1] <- m.1
                tmp[, 2] <- m.2
              }
              if (is.element(receptors, tmp[, 2])) {
                siz <- c(siz, rowMeans(data[tmp[tmp[,2] %in% receptors, 1], 
                                            cluster == as.numeric(which(c.names %in% colnames(tmp)[1])),
                                            drop = FALSE] # advised by duanmuseventeen
                                       )
                         )
                nam.tmp <- rep(colnames(tmp)[1], nrow(tmp[tmp[, 
                                                              2] %in% receptors, ]))
                colnames(tmp)[seq_len(2)] <- c("a.gn", 
                                               "b.gn")
                tmp[, 1] <- paste(nam.tmp[1], tmp[, 
                                                  1], sep = "-")
                add.net <- rbind(add.net, tmp[tmp[, 
                                                  2] %in% receptors, ])
                nam <- c(nam, nam.tmp)
              }
            }
          }
          if (is.null(add.net) == FALSE) {
            add.net <- cbind(add.net[, seq_len(2)], 
                             nam, location = "extra", type = "control")
          }
        }
        net.tmp <- cbind(net.n[, seq_len(2)], nam = rep(coi, 
                                                        nrow(net.n)), location = "intra", type = net.n$type)
        net.f <- rbind(add.net, net.tmp)
        if (species == "mus musculus") {
          net.f[, 1] <- Hs2mm[net.f[, 1]]
          net.f[, 2] <- Hs2mm[net.f[, 2]]
        }
        g.net <- graph_from_data_frame(net.f, directed = TRUE)
        g.net.tmp <- as.undirected(g.net)
        y <- shortest_paths(g.net.tmp, receptors, V(g.net.tmp))
        g.net.tmp <- graph_from_data_frame(net.f[, seq_len(2)], 
                                           directed = FALSE)
        V(g.net.tmp)$status <- "pw.related"
        V(g.net.tmp)$status[which(unique(c(net.f$a.gn, 
                                           net.f$b.gn)) %in% receptors)] = "gene.of.interest"
        V(g.net.tmp)$status[which(unique(c(net.f$a.gn, 
                                           net.f$b.gn)) %in% net.f$a.gn[net.f$location == 
                                                                          "extra"])] = "ligand"
        E(g.net.tmp)$int.type <- as.character(net.f$type)
        y <- unlist(lapply(y$vpath, function(x) length(x)))
        y <- y - 1
        if (sum(y == -1) > 0 & connected == FALSE) {
          y[y == -1] <- 1
        }
        if (sum(y == -1) > 0 & connected == TRUE) {
          nam.tmp <- unique(c(net.f$a.gn, net.f$b.gn))[y == 
                                                         -1]
          net.f <- net.f[!net.f$a.gn %in% nam.tmp & 
                           !net.f$b.gn %in% nam.tmp, ]
          y <- y[y != -1]
          g.net <- graph_from_data_frame(net.f, directed = TRUE)
        }
        names(y) <- unique(c(net.f$a.gn, net.f$b.gn))
        y[net.f$a.gn[net.f$location == "extra"]] = -1
        x = vector("numeric", length = length(y))
        names(x) <- unique(c(net.f$a.gn, net.f$b.gn))
        if (sum(net.f$location == "extra")) {
          x[net.f$a.gn[net.f$location == "extra"]] <- (seq(0, 
                                                           length(net.f$a.gn[net.f$location == "extra"]) * 
                                                             4 - 1, 4) - max(seq(0, length(net.f$a.gn[net.f$location == 
                                                                                                        "extra"]) * 4 - 1, 4))/2)
        }
        x[receptors] <- 0
        for (j in seq_len(max(y))) {
          if (sum(y == j) > 1) {
            x[y == j] <- seq(-10, 10, 19/sum(y == j))[seq_len(sum(y == 
                                                                    j))]
          }
        }
        l <- cbind(x, -y)
        for (i in seq_len(max(y))) {
          if (sum(y == i) > 1) {
            index <- which(y == i)
            l[index[seq(1, length(index), 2)], 2] <- l[index[seq(1, 
                                                                 length(index), 2)], 2] - 0.2
            l[index[seq(2, length(index), 2)], 2] <- l[index[seq(2, 
                                                                 length(index), 2)], 2] + 0.2
          }
        }
        rownames(l) <- unique(c(net.f$a.gn, net.f$b.gn))
        if (sum(net.f$location %in% "intra") == 0) {
          V(g.net)$size[y == -1] <- (log(siz) + 5) * 
            4
          V(g.net)$size[V(g.net)$size < exp(-5)] <- exp(-5)
        } else {
          V(g.net)$size <- (log(c(siz, rowMeans(data.tmp[unique(c(net.n$a.gn, 
                                                                  net.n$b.gn)), ]))) + 5) * 4
          V(g.net)$size[V(g.net)$size < exp(-5)] = exp(-5)
        }
        V(g.net)$size[which(unique(c(net.f$a.gn, net.f$b.gn)) %in% 
                              receptors)] <- 35
        V(g.net)$shape <- c("circle")
        V(g.net)$shape[which(unique(c(net.f$a.gn, net.f$b.gn)) %in% 
                               receptors)] <- c("rectangle")
        V(g.net)$color <- c("lightcyan")
        V(g.net)$color[which(unique(c(net.f$a.gn, net.f$b.gn)) %in% 
                               receptors)] <- c("indianred1")
        V(g.net)$color[which(unique(c(net.f$a.gn, net.f$b.gn)) %in% 
                               net.f$a.gn[net.f$location == "extra"])] <- c("lightyellow1")
        V(g.net)$label.color <- "black"
        V(g.net)$label.dist <- 0 # advised by duanmuseventeen
        V(g.net)$label.dist[which(!unique(c(net.f$a.gn, net.f$b.gn)) %in% receptors)] <- 1
        for (i in unique(l[, 2])) {
          if (i != 0) {
            V(g.net)$label.degree[l[, 2] == i] <- rep(c(pi/2, 
                                                        -pi/2), length(y))[seq_len(sum(l[, 2] == 
                                                                                         i))]
          }
        }
        V(g.net)$label.dist[which(unique(c(net.f$a.gn, 
                                           net.f$b.gn)) %in% receptors)] <- 0
        V(g.net)$label.degree[which(unique(c(net.f$a.gn, 
                                             net.f$b.gn)) %in% receptors)] <- 0
        E(g.net)$arrow.mode <- as.numeric(grepl("control", 
                                                net.f$type)) * 2
        E(g.net)$arrow.size <- rep(0.4, nrow(net.f))
        E(g.net)$color <- "gray"
        g.pw.names <- strsplit(net.n$pathway, ";")
        g.pw.names <- g.pw.names[g.pw.names != "added"]
        g.pw.sizes <- table(unlist(g.pw.names))
        N <- nrow(PwC_ReactomeKEGG)
        n <- nrow(net.n)
        pw.table <- foreach(pw = names(g.pw.sizes), 
                            .combine = rbind) %do% {
                              K <- pw.sizes[pw]
                              if (K <= max.pw.size) {
                                k <- g.pw.sizes[pw]
                                pval <- 1 - phyper(q = k - 1, m = K, n = N - 
                                                     K, k = n)
                                data.frame(pathway = pw, in.pw = k, pw.size = K, 
                                           pval = pval, stringsAsFactors = FALSE)
                              }
                              else NULL
                            }
        par(mfrow = c(2, 1))
        par(las = 2)
        par(mar = c(0, 0, 1, 0))
        plot(g.net, layout = l, main = coi)
        if (is.null(pw.table) == FALSE) {
          rawp <- pw.table$pval
          if (length(rawp) == 1) {
            adj <- rawp
            qval <- adj
          }
          else {
            adj <- mt.rawp2adjp(rawp, "BH")
            qval <- adj$adjp[order(adj$index), "BH"]
          }
          pw.table <- cbind(pw.table, data.frame(qval = qval))
          rownames(pw.table) <- NULL
          pw.table <- pw.table[pw.table$qval < 0.05, 
          ]
          if (nrow(pw.table) > 0) {
            pw.table <- pw.table[order(pw.table$in.pw, 
                                       decreasing = FALSE), ]
            nc <- max(nchar(as.character(pw.table$pathway))) * 
              0.3
            mtext(paste(receptors, "related pathways"), 
                  side = 1, line = 2, las = FALSE)
            par(mar = c(2, nc, 4, 2))
            barplot((pw.table$in.pw), horiz = TRUE, 
                    names.arg = (paste(pw.table$pathway, "*")), 
                    cex.names = 0.7, col = "azure2", border = "gray60")
          }
          else {
            if (verbose == TRUE) {
              cat("No significant associated pathway", 
                  fill = TRUE)
            }
          }
        }
        else {
          if (verbose == TRUE) {
            cat(paste("No associated genes downstream", 
                      receptors, "in", coi), fill = TRUE)
          }
        }
        if (write == TRUE) {
          if (is.null(pw.table) == FALSE) {
            write.table(pw.table[order(pw.table$qval), 
            ], file = paste0("./networks/intracell_network_pathway_analysis_", 
                             coi, "-", receptors, ".txt"), sep = "\t", 
            quote = FALSE, row.names = FALSE)
          }
          write.graph(g.net.tmp, paste0("./networks/intracell_network_", 
                                        coi, "-", receptors, ".graphml"), format = "graphml")
        }
        res[[qq]] <- net.f
      }
      else {
        cat("No interactions found.", fill = TRUE)
      }
    }
  }
  par(opar)
  return(res)
}
