# sankey with ggsankey
# https://github.com/davidsjoberg/ggsankey
# https://mp.weixin.qq.com/s/uN3pypOU1F70l86y8iujKA

require(dplyr)
require(clusterProfiler)
require(ggplot2)
require(networkD3)
require(AnnotationDbi)
require(ggsankey)

# Load data---------------------------------------------------------------------
load("35.pathway.Rdata")

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

dat.sankey.ml <- dat.sankey %>% 
  make_long(pathway, symbol)

library(RColorBrewer)

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

ggplot(dat.sankey.ml, aes(
  x = x, 
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  label = node,
  value = 1)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = "#eeffff55", 
              node.fill = mycol,
              smooth = 8,
              width = 0.08,
              type = "alluvial") +
  geom_sankey_text(size = 3.2, color = "black",type = "alluvial")+
  theme_void() +
  theme(legend.position = 'none')