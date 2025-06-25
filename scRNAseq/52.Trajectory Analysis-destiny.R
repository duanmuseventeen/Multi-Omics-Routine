library(Seurat)
library(destiny)

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

# local-------------------------------------------------------------------------
# https://www.bioconductor.org/packages/release/bioc/vignettes/destiny/inst/doc/Diffusion-Maps.html

# The data necessary to create a diffusion map with our package is a a cell×gene matrix or data.frame, or alternatively an ExpressionSet (which has a gene×cell exprs matrix). In order to create a DiffusionMap object, you just need to supply one of those formats as first parameter to the DiffusionMap function. In the case of a data.frame, each floating point column is interpreted as expression levels, and columns of different type (e.g. factor, character or integer) are assumed to be annotations and ignored. Note that single-cell RNA-seq count data should be transformed using a variance-stabilizing transformation (e.g. log or rlog); the Ct scale for qPCR data is logarithmic already (an increase in 1 Ct corresponds to a doubling of transcripts).
# In order to create a diffusion map to plot, you have to call DiffusionMap, optionally with parameters. If the number of cells is small enough (< ~1000), you do not need to specify approximations like k (for k nearest neighbors).
# If you started reading here, execute data(guo.norm) to load the dataset that was created in the previous section.

# dm <- DiffusionMap(as.SingleCellExperiment(fig, assay = "RNA"), verbose = TRUE)
dm <- DiffusionMap(t(as.matrix(fig@assays$RNA$data)), verbose = TRUE)
plot(dm, pch = 20)
plot(dm, col = factor(fig$cell_type, 
                      levels = c("A","B","C","D"), 
                      labels = colors_list[c(1,3,5,7)]) %>% as.character, 
     pch = 20)

plot(dm, 1:2, legend_main = 'Cell stage')

library(ggplot2)
# The vectors of DC are stored in dm@eigenvectors
ggplot(dm, aes(DC1, DC2, colour = factor(num_cells))) +
  geom_point() + 
  scale_color_cube_helix()
# global------------------------------------------------------------------------
sigmas <- find_sigmas(t(as.matrix(fig@assays$RNA$data)), verbose = TRUE)
par(lwd = 3)
plot(sigmas,
     col           = palette()[[1]],
     col_highlight = palette()[[4]],
     col_line      = palette()[[6]])

dm_global <- DiffusionMap(
  as.SingleCellExperiment(fig, assay = "RNA"), 
  sigmas, 
  verbose = FALSE)
plot(dm_global, 
     col = factor(fig$cell_type, 
                  levels = c("Neu_S100A12","Neu_NLRP3","Neu_CCL4","Neu_IFIT1"),
                  labels = colors_list[c(1,3,5,7)]) %>% as.character, 
     pch = 20)