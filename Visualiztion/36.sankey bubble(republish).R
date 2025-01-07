# sankey with ggsankey
# https://github.com/davidsjoberg/ggsankey
# https://mp.weixin.qq.com/s/uN3pypOU1F70l86y8iujKA

require(dplyr)
require(clusterProfiler)
require(ggplot2)
require(networkD3)
require(AnnotationDbi)
require(ggsankey)
require(AnnotationDbi)
require(RColorBrewer)
# Load data---------------------------------------------------------------------
dat <- readxl::read_excel("C:/D/Labmates/金思佳/2025-01-06 sankey/GSE183464_NEU1 for R.xlsx",
                          sheet = "Sheet4")

# ORA---------------------------------------------------------------------------
geneList <- bitr(dat$gene_name, fromType = "SYMBOL", toType = "ENTREZID", "org.Hs.eg.db")

# ERA_KEGG <- enrichKEGG(
#   gene = geneList$ENTREZID,
#   organism = "hsa",
#   keyType = "kegg",
#   pvalueCutoff = 1, 
#   qvalueCutoff = 1)
# ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

ERA_GO <- enrichGO(
  gene = geneList$ENTREZID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID", 
  ont = "ALL",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = T)
ERA_GO <- setReadable(ERA_GO, 'org.Hs.eg.db', 'ENTREZID')
# graphics----------------------------------------------------------------------
# pathway <- ERA_GO@result %>% 
#   filter(Description %in% c("glycosylation",
#                             "glycolipid metabolic process",
#                             "membrane lipid metabolic process",
#                             "protein O-linked glycosylation",
#                             "sialylation",
#                             "N-glycan processing",
#                             "monoatomic anion transport",
#                             "glycoprotein catabolic process",
#                             "lipid glycosylation",
#                             "N-acetylglucosamine metabolic process"
#   ))

symbol <- str_split(pathway$geneID, "/")
names(symbol) <- pathway$Description

level.p <- symbol %>% unlist2 %>% names %>% unique %>% sort
level.s <- symbol %>% unlist %>%  unique %>% sort

dat.sankey <- data.frame(
  pathway = symbol %>% unlist2 %>% names,
  symbol  = symbol %>% unlist,
  value   = 1,
  stringsAsFactors = FALSE,
  row.names = NULL
) %>%
  mutate(pathway = factor(pathway, levels = level.p),
         symbol = factor(symbol, levels = level.s))

dat <- dat %>% 
  mutate(Order = as.character(Order)) %>% 
  mutate(Order = case_when(
    nchar(Order) == 1 ~ paste0("00", Order),
    nchar(Order) == 2 ~ paste0("0", Order),
    nchar(Order) == 3 ~ paste0("", Order)
    )) %>% 
  mutate(symbol = gene_name,
         order_sybmol = paste0(Order,"::",gene_name))

dat.sankey <- dat.sankey %>% 
  left_join(dat %>% dplyr::select(symbol, order_sybmol), by = "symbol") %>% 
  transmute(pathway = pathway, symbol = order_sybmol)

dat.sankey.ml <- dat.sankey %>% 
  make_long(pathway, symbol)

mycol <- c(
  # scales::muted(rainbow(12)),
  rainbow(12),
  colorRampPalette(c("#F37252", "#FFF0DC", "#779FC6"))(48)
  )

# the length of values are 2-fold of nodes
# values ranks from left-bottom and the matched node in the right, then,
# left-bottom second and the matched node, repeat
# in this example, set value = 1 makes all symbols have the same width

# fun first rank the nodes, then colors are assigned

p <- ggplot(dat.sankey.ml, aes(
  x = x, 
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  label = node,
  value = 1)) +
  geom_sankey(flow.alpha = 0.5,
              # flow.fill = 'grey',
              # flow.color = "#eeffff55", 
              # node.fill = mycol,
              smooth = 8,
              width = 0.08,
              type = "alluvial") +
  geom_sankey_text(size = 3.2, color = "black",type = "alluvial")+
  theme_void() +
  theme(legend.position = 'none')

ggsave(p, filename = "sankey.pdf", width=10, height=10, units="in")
