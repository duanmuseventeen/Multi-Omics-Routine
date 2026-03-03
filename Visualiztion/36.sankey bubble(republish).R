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
dat <- readxl::read_excel("dat for R.xlsx")

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
pathway <- ERA_GO@result %>% 
  filter(Description %in% c("glycosylation",
                            "glycolipid metabolic process",
                            "membrane lipid metabolic process",
                            "protein O-linked glycosylation",
                            "sialylation",
                            "N-glycan processing",
                            "monoatomic anion transport",
                            "glycoprotein catabolic process",
                            "lipid glycosylation",
                            "N-acetylglucosamine metabolic process"
  ))

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

# 2026-03-03 example2 ----
require(dplyr)
require(ggplot2)
require(AnnotationDbi)
require(ggsankey)
require(AnnotationDbi)
require(RColorBrewer)
  
  head(node)
  # gene        coef                                 id         db subtype symbol
  # 1 LPAR6 -0.12065943 HALLMARK_INTERFERON_ALPHA_RESPONSE Hallmark50          LPAR6
  # 2 ERCC1 -0.11184625                         GO:0002449         GO      BP  ERCC1
  # 3 ERCC1 -0.11184625                         GO:0002460         GO      BP  ERCC1
  # 4 NFIL3 -0.09102736   HALLMARK_TNFA_SIGNALING_VIA_NFKB Hallmark50          NFIL3
  # 5  BTG1 -0.08676659 HALLMARK_INTERFERON_GAMMA_RESPONSE Hallmark50           BTG1
  # 6  BTG1 -0.08676659   HALLMARK_TNFA_SIGNALING_VIA_NFKB Hallmark50           BTG1
  
  head(output)
  # gene        coef
  # LPAR6 LPAR6 -0.12065943
  # ERCC1 ERCC1 -0.11184625
  # NFIL3 NFIL3 -0.09102736
  # BTG1   BTG1 -0.08676659
  # EXT1   EXT1 -0.08488127
  # YAP1   YAP1 -0.08480225
  
  dat.sankey <- node %>% arrange(desc(coef))
  dat.sankey$gene <- factor(dat.sankey$gene, levels = unique(dat.sankey$gene))
  dat.sankey <- dat.sankey %>% dplyr::select(gene, id, db)
  
  dat <- output %>% 
    arrange(coef) %>% 
    mutate(Order = as.character(row_number())) %>% 
    mutate(Order = case_when(
      nchar(Order) == 1 ~ paste0("0", Order),
      nchar(Order) == 2 ~ paste0("", Order)
    )) %>% 
    mutate(symbol = gene,
           order_sybmol = paste0(Order,"::",gene))
  
  dat.sankey <- dat.sankey %>% 
    mutate(value = 1, symbol = gene, pathway = id) %>% 
    left_join(dat %>% mutate(symbol = gene) %>% dplyr::select(symbol, order_sybmol), by = "symbol") %>% 
    transmute(pathway = pathway, symbol = order_sybmol, value = value)
  
  dat.sankey.ml <- dat.sankey %>% 
    make_long(pathway, symbol)
  
  values <- dat.sankey %>% 
    group_by(symbol) %>%
    mutate(symbol_value = 1/n()) %>%
    mutate(pathway_value = symbol_value) %>% 
    ungroup() %>% 
    arrange(pathway) %>% 
    arrange(desc(symbol)) %>% 
    dplyr::select(symbol_value, pathway_value) %>% 
    as.matrix %>% t %>% 
    as.list %>% unlist
  
  mycol <- c(
    # scales::muted(rainbow(12)),
    rep('#FFE6B7', 31), # GO
    rep('#C7B8BD', 6), # HALL
    rep('#AADCE0', 2), # REACTOME 
    rep('#f4a494', 7), # KEGG
    colorRampPalette(c('#1E466E', "white","red3" ))(71)
  )
  
  ggplot(dat.sankey.ml, aes(
    x = x, 
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = node,
    value = values)) +
    geom_sankey(flow.alpha = 0.5,
                # flow.fill = 'grey',
                # flow.color = "#eeffff55", 
                node.fill = mycol,
                smooth = 8,
                width = 0.08,
                type = "alluvial") +
    geom_sankey_text(size = 3.2, color = "black",type = "alluvial")+
    theme_void() +
    theme(legend.position = 'none')
