# Subcellular localization
# download from UNIPROT---------------------------------------------------------
# https://www.uniprot.org/
annot_cc <- data.table::fread("UNIPROT_Species.tsv")

colnames(gene2annot)[2] <- "HSA Gene Symbol" # a df obj including genes 2 annotate

annot_cc
[1] "Entry"	"Entry name"	"Status"	"Protein names"	"Gene names"	"Organism"	"Subcellular location [CC]"
colnames(annot_cc)[1] <- "HSA Gene Symbol" # Entry

gene2annot %>%
  dplyr::left_join(annot_cc, by = "HSA Gene Symbol") %>% 
  export::table2excel("genes with Subcellular localization.xlsx")
# download from ProteinAtlas (only human) --------------------------------------
# https://www.proteinatlas.org/humanproteome/subcellular/data#locations
annot_cc <- data.table::fread("subcellular_location.tsv")

colnames(gene2annot)[2] <- "HSA Gene Symbol" # a df obj including genes 2 annotate
colnames(annot_cc)[2] <- "HSA Gene Symbol"

gene2annot %>%
  dplyr::left_join(annot_cc, by = "HSA Gene Symbol") %>% 
  export::table2excel("genes with Subcellular localization.xlsx")