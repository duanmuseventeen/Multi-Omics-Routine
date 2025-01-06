# sankey with networkD3
# https://r-graph-gallery.com/323-sankey-diagram-with-the-networkd3-library.html

require(dplyr)
require(clusterProfiler)
require(ggplot2)
require(networkD3)
require(AnnotationDbi)
require(stringr)

# Load data---------------------------------------------------------------------
load("35.pathway.Rdata")
# prepare-----------------------------------------------------------------------
# symbol <- str_split(pathway$geneID, "/",simplify = TRUE)
symbol <- str_split(pathway$geneID, "/")
names(symbol) <- pathway$Description

dat.sankey <- data.frame(
  pathway = symbol %>% unlist2 %>% names,
  symbol  = symbol %>% unlist,
  value   = 1,
  stringsAsFactors = FALSE,
  row.names = NULL
) %>% 
  group_by(symbol) %>% 
  mutate(value = value / n())
annot <- data.frame(
  name = unique(c(dat.sankey$pathway, dat.sankey$symbol))
) %>% 
  mutate(num = name %>% factor(levels = name) %>% as.numeric) %>% 
  mutate(num = num - 1)

dat.sankey.num <- dat.sankey %>% 
  mutate(name = pathway) %>% 
  left_join(annot, by = "name") %>% 
  transmute(pathway = pathway, symbol = symbol , source = num, value = value) %>% 
  mutate(name = symbol) %>% 
  left_join(annot, by = "name") %>% 
  transmute(source = source, target = num, value = value)

# graphics----------------------------------------------------------------------
# Thus we can plot it
p <- networkD3::sankeyNetwork(
  Links = dat.sankey.num, 
  Nodes = annot %>% dplyr::select(name), 
  Source = "source",
  Target = "target", 
  Value = "value",
  NodeID = "name",
  # units = "TWh",
  # sinksRight=FALSE, 
  fontSize = 12, 
  nodeWidth = 30)

# save--------------------------------------------------------------------------
library(htmlwidgets)
saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyEnergy.html"))

saveNetwork(p, "F:/sankey_diagram.html")
  