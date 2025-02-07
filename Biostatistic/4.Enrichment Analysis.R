# https://blog.csdn.net/qq_50898257/article/details/120588222


# N: 基因组所有基因、所有分析的基因
# x: 差异表达基因集中有功能F的基因
# M: N中具有某种功能（F）的基因总数
# K: 差异表达基因

enrich <- function(N,x,M,K){
  res <- 0
  for (i in 0:(x - 1)) {
    res <- res + choose(M,i) * choose((N - M),(K - i)) / choose(N, K)
  }
  p.value <- 1 - res
  return(p.value)
}

# example 1 --------------------------------------------------------------------
# TermID	Term	bg_pro_num	bg_term_num	fg_pro_num	fg_term_num	Ratio	pvalue	Enrichment	FDR	IDs	GeneSymbol	Link
# GO:0043066	negative regulation of apoptotic process	88	15	29	11	0.733333333	0.000553211	3.257109327	0.392779689	GHRL_MOUSE,ERBB4_MOUSE,FOXO1_MOUSE,TNR6_MOUSE,TGFA_MOUSE,GLUC_MOUSE,EPCAM_MOUSE,SCF_MOUSE,TNFA_MOUSE,PRDX5_MOUSE,IL6_MOUSE	Ghrl,Erbb4,Foxo1,Fas,Tgfa,Gcg,Epcam,Kitlg,Tnf,Prdx5,Il6	https://www.ebi.ac.uk/QuickGO/term/GO:0043066
enrich(N = 88, x = 11, M = 15, K = 29)
[1] 0.0005532108

phyper(q = 11 - 1, m = 29, n = 88 - 29, k = 29, lower.tail = FALSE)
# [1] 0.322076
# example 2 --------------------------------------------------------------------
# ONTOLOGY         ID                                              Description
# GO:0050731       BP GO:0050731 positive regulation of peptidyl-tyrosine phosphorylation
# GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# GO:0050731     16/87 182/28905 2.504858e-19 7.183933e-16 3.193035e-16
# geneID
# GO:0050731 Il6/Tgfa/Tnf/Gfra1/Erbb4/Gdnf/Csf2/Ntf3/Pdgfb/Yes1/Tgfb1/Ccl5/Epo/Cntn1/Hgf/Il5
# Count
# GO:0050731    16
enrich(28905, 16, 182, 87)

phyper(q = 16 - 1, m = 182, n = 28905 - 182, k = 87, lower.tail = FALSE)
# [1] 2.504858e-19

# category                         subcategory
# mmu04060 Environmental Information Processing Signaling molecules and interaction
# ID                                                         Description
# mmu04060 mmu04060 Cytokine-cytokine receptor interaction - Mus musculus (house mouse)
# GeneRatio  BgRatio      pvalue     p.adjust       qvalue
# mmu04060     24/66 294/9574 3.80724e-20 7.119539e-18 4.368307e-18
# geneID
# mmu04060 Il6/Cxcl1/Il1b/Il1a/Ccl2/Tnf/Tnfrsf12a/Ccl3/Cxcl9/Fas/Tnfrsf11b/Csf2/Acvrl1/Il23r/Il10/Ccl20/Il17f/Il17a/Tgfb1/Ccl5/Epo/Il5/Eda2r/Tnfsf12
# Count
# mmu04060    24
enrich(N = 9574, x = 24, M = 294, K = 66)

phyper(q = 24 - 1, m = 294, n = 9574 - 294, k = 66, lower.tail = FALSE)
# [1] 3.80724e-20

if(FALSE){
  # fisher.test(
  #   matrix(c( 66 - 24, 24, 9574 - 294, 294),nrow = 2)
  #   )$p.value
  # # [1] 2.078765e-19
  # fisher.test(
  #   matrix(c( 24, 66 - 24, 294, 9574 - 294),nrow = 2)
  # )$p.value
  # # [1] 2.078765e-19
  # fisher.test(
  #   matrix(c( 294, 9574 - 294, 24, 66 - 24),nrow = 2)
  # )$p.value
  # # [1] 2.078765e-19
  # 
  # chisq.test(
  #   matrix(c( 294, 9574 - 294, 24, 66 - 24),nrow = 2)
  # )$p.value
  # [1] 3.265787e-49
}

