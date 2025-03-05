# Predict miRNA-mRNA
# Please download miRNA-mRNA files from databases manually

# Database======================================================================
# 【TargetScanHuman】https://www.targetscan.org/vert_72/
# 【RNA22 v2 microRNA target detection】https://cm.jefferson.edu/rna22/Interactive
#  https://cm.jefferson.edu/rna22-full-sets-of-predictions/
#  https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid
# 【miRcode】http://www.mircode.org/index.php
# 【miRDB】https://mirdb.org/cgi-bin/search.cgi
# 【Unavailable】https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php
# 【DIANA】http://www.microrna.org/microrna/home.do
#  https://pictar.mdc-berlin.de/
# 【Page not found】https://genie.weizmann.ac.il/pubs/mir07/mir07_dyn_data.html
# 【plant】https://www.zhaolab.org/psRNATarget/
# 【plant】http://omicslab.genetics.ac.cn/psRobot/target_prediction_1.php
#*【Starbase: API】https://rnasysu.com/encori/
#* https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/php/index.php
#*【404】http://mirtarbase.mbc.nctu.edu.tw/php/index.php


targetmiRNA <- c("")

# starBase [API]---------------------------------------------------------------
# https://rnasysu.com/encori/tutorialAPI.php
dat <- readxl::read_excel("starBase.xlsx") # Result from starBase

match.mRNA.starBase <- dat %>%
  filter(miRNAname == "hsa-miR-221-3p") %>%
  dplyr::select(geneName) %>%
  unlist %>%
  unique

# mircode----------------------------------------------------------------------
mircode <- data.table::fread("mircode_highconsfamilies.txt.gz")

mircode$microrna[which(str_detect(mircode$microrna, "132"))] %>% unique
## miRTarBase-------------------------------------------------------------------
miRTarBase <- readxl::read_excel("miRTarBase_MTI.xlsx") %>% 
  filter(`Species (miRNA)` == "Homo sapiens")
match.mRNA.miRTarBase <- miRTarBase %>% 
  filter(miRNA %in% targetmiRNA) 
## TarBase-v9----------------------------------------------------------------------
TarBase <- data.table::fread("Homo_sapiens_TarBase-v9.tsv.gz")

match.mRNA.TarBase <- TarBase %>%
  filter(species == "Homo sapiens") %>%
  filter(mirna_name %in% targetmiRNA) %>%
  dplyr::select(gene_name) %>%
  unlist %>%
  unique
## targetscan----------------------------------------------------------------------
targetscan <- data.table::fread("Predicted_Targets_Context_Scores.default_predictions.txt")
match.mRNA.targetscan <- targetscan %>%
  filter(miRNA %in% targetmiRNA) %>%
  dplyr::select(`Gene Symbol`) %>%
  unlist %>%
  unique