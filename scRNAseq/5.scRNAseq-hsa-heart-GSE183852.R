# Single-cell transcriptomics reveals cell-type-specific diversification in human heart failure----
# https://pubmed.ncbi.nlm.nih.gov/35959412/
# GSE183852
  
# RefMerge@meta.data$Names is cell type
# RefMerge@meta.data$condition is group information

# Load data
# Download from https://explore.data.humancellatlas.org/projects/135f7f5c-4a85-4bcf-9f7c-4f035ff1e428
load(GSE183852_DCM_Integrated.Robj)

# Visualization
DimPlot(RefMerge, reduction = "harmony",
              group.by = "Names", 
              label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(RefMerge, features = c(Targetgene))
VlnPlot(object = RefMerge, features = Targetgene, slot = "counts", log = TRUE, group.by = "Names")

# Subset
RefMerge <- RegroupIdents(RefMerge, metadata = "Names")
CM <- subset(RefMerge, subset = Names == "Cardiomyocytes")

# Differential Analysis
Targetgene <- CM@assays$RNA$counts[rownames(CM@assays$RNA$counts) == Targetgene]
group <- CM@meta.data$condition

mean(Targetgene[group == "Donor"])
mean(Targetgene[group != "Donor"])
wilcox.test(Targetgene ~ group)

dat2pseudo <- data.frame(
    ID = CM@meta.data$orig.ident,
    disease = CM@meta.data$condition,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>% 
    filter(!duplicated(ID)) %>% 
    mutate(Targetgene = NA)
  
  for (i in unique(CM@meta.data$orig.ident)) {
    tmp <- subset(x = CM, subset = orig.ident == i) 
    dat2pseudo$Targetgene[dat2pseudo$ID == i] <- 
      tmp@assays$RNA$counts[rownames(tmp@assays$RNA$counts) == Targetgene] %>% sum
    print(i)
    }
  
  dat2pseudo %>% group_by(disease) %>% mutate(mean = mean(Targetgene)) %>% distinct(mean)
  wilcox.test(Targetgene ~ disease, data = dat2pseudo)
