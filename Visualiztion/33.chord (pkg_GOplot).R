# chord plot
rm(list = ls());gc()

# Usually, this plot are performed from a sequencing data processed by limma, DESeq2
# or edgeR. Here, I just use GOChord to perform a chord plot to display the association
# between genes and pathways. Of course, you can use this way to display any
# relationship between category and units, for example metabolites.

# load pkgs---------------------------------------------------------------------
require(ggtext)
require(dplyr)
require(clusterProfiler)
require(enrichplot)
require(GOplot)
require(stringr)

# graph-------------------------------------------------------------------------
# I create the pathway and genelist files according to tutorial. For our purpose
# some information are not necessary, they are left empty
pathway <- readxl::read_excel("pathway.xlsx")
genelist <- readxl::read_excel("genelist.xlsx")

circ <- circle_dat(pathway, genelist)

chord <- chord_dat(data = circ, genes = genelist %>% mutate(ID = toupper(ID)) %>% as.data.frame)

pdf("chord.pdf", width = 9, height = 10)
GOChord(chord, 
        space = 0.02, 
        gene.order = 'logFC', 
        gene.space = 0.25, 
        gene.size = 5,
        ribbon.col = c("#F37252", "#F7935A", "#FCB461", "#FED297", "#FFF0DC",
                       "#E3F5F6", "#ABE0E4", "#7FC6D3", "#779FC6"))
dev.off()
