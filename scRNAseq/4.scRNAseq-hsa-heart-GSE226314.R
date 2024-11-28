# Defining cardiac functional recovery in end-stage heart failure at single-cell resolution.----
# https://pubmed.ncbi.nlm.nih.gov/37583573/
# https://explore.data.humancellatlas.org/projects/c0fecf0b-af86-41b8-ba82-d5fd81b7542a
# GSE226314

# Sample selection clinical phenotyping
# The recovery cohort of patients was selected to match for the following variables to the best extent: ejection fraction at the time of LVAD implant, sex, age and clinical risk profile. Age-matched and sex-matched donors were then pulled from the Washington University in St. Louis School of Medicine biobank repository. Within the LVAD cohort, patients were assigned as ‘reverse remodeled’ or ‘unloaded’. Donors were selected to age and sex match with the LVAD samples. RR and U samples were chosen such that ejection fraction at the time of LVAD implantation was not different between the two groups.

Supplementary Table 1: Demographic Details of samples							

Sample	Condition	Age	Sex	BMI	EF	HTN	DM
H_ZC-11-292	Donor	68	male	24.7		1	0
TWCM-11-41	Donor	60	male	23.3	62.5	1	1
TWCM-11-74	Donor	65	male	23.7	55	1	0
TWCM-11-78	Donor	27	male	24.54	47.5	1	1
TWCM-11-82	Donor	44	male	25.81	60	0	0
TWCM-11-103	Donor	46	male	29.9	60	1	0
TWCM-11-192	Donor	63	male	25.6	70	0	0
TWCM-13-1	Donor	61	female	42.76	65	1	1
TWCM-13-80	Donor	46	male	32.87	65	0	0
TWCM-13-104	Donor	48	female	22.88	50	0	0
TWCM-13-152	Donor	35	male	23.81	74	0	0
TWCM-13-168	Donor	66	male	26.5	66	1	1
TWCM-13-192	Donor	21	male	30.56	65	0	1
TWCM-14-173	Donor	63	female	22.51	60	0	0
DONOR AVERAGE		50.92857143		27.10285714	61.53846154		

190	Responder	18	Male	22.18	12	0	0
229	Responder	60	Male	26.904	35	1	0
239	Responder	44	Male	29.758	15	0	0
296	Responder	60	Female	26.968	25	1	1
359	Non Responder	46	Male	25.448	15	0	0
363	Non Responder	64	Male	26.646	24	1	0
370	Non Responder	67	Female	24.242	13	0	0
371	Non Responder	61	Male	24.205	11	1	1
373	Non Responder	69	Male	23.755	10	0	0
376	Non Responder	37	Male	22.498	10	0	0
397	Non Responder	47	Female	26.523	15	1	0
410	Non Responder	58	Male	27.335	23	0	0
463	Responder	29	Male	38.729	5	1	0
BASELINE LVAD AVERAGE		50.76923077		26.55315385	16		

# Load data
sce.all <- readRDS("GSE226314_global.rds")

# Subset data
# Because this analysis aims to compare the RNA expression between normal and failing
# hearts. The sample after implantation are filtered.

sce.all <- RegroupIdents(sce.all, metadata = "condition")
# Don't use idents para to filter, even it is used in tutorial.
sce.all <- subset(x = sce.all, subset = condition %in% c("Donor", "NRpre", "Rpre"))

# Visualization
DimPlot(sce.all, reduction = "umap",
        group.by = "cell.type",
        label = TRUE, pt.size = 0.5)
DimPlot(sce.all, reduction = "umap",
        group.by = "condition",
        label = TRUE, pt.size = 0.5)
FeaturePlot(sce.all, features = c(Targetgene), raster = FALSE)
VlnPlot(object = sce.all, features = Targetgene, slot = "counts", log = TRUE, group.by = "cell.type")

# Subset Caridal Myocytes
sce.all <- RegroupIdents(sce.all, metadata = "cell.type")
CM <- subset(x = sce.all, subset = cell.type == "Cardiomyocyte") 

CM@meta.data$disease <- CM@meta.data$condition
CM@meta.data$disease[CM@meta.data$disease == "Donor"] <- "Normal"
CM@meta.data$disease[CM@meta.data$disease != "Normal"] <- "Heart Failure"

VlnPlot(object = CM, features = Targetgene, slot = "counts", log = TRUE, group.by = "disease")

# Donor <- subset(CM, subset = condition %in% c("Donor"))
# HF <- subset(CM, subset = condition %in% c("NRpre", "Rpre"))
Donor <- subset(CM, subset = disease == "Normal")
HF <-    subset(CM, subset = disease != "Normal")

# Differential Analysis
Donor@assays$RNA$counts[rownames(Donor@assays$RNA$counts) == Targetgene,] %>% unlist %>% mean
HF@assays$RNA$counts[rownames(HF@assays$RNA$counts) == Targetgene,] %>% unlist %>% mean
wilcox.test(Donor@assays$RNA$counts[rownames(Donor@assays$RNA$counts) == Targetgene,] %>% unlist,
            HF@assays$RNA$counts[rownames(HF@assays$RNA$counts) == Targetgene,] %>% unlist)