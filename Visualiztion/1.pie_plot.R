require(dplyr)
require(ggplot2)

df <- data.frame(
  x = c("yaocai","yaocai","chengyao","chengyao") %>% factor(levels = c("yaocai","chengyao")),
  y = c(80 / (80 + 536), 536 / (80 + 536), 558 / (558+1049), 1049/(558+1049)),
  col = c("a","b","c","d"),
  # y = c(),
  labelx = c(80, 536, 558, 1049) %>% as.factor,
  stringsAsFactors = F
)

p <- ggplot(data = df, aes(x = x, y = y, fill = col)) + 
  geom_col(col = "white") + 
  # ggrepel::geom_text_repel() +
  scale_fill_manual(values = c("#ECA8A9","#ddeedd","#F7C97E","#ddeeff")) +
  guides(fill = "none", col = "none") +
  scale_x_discrete(labels = c('','')) +
  # theme(aspect.ratio=1) + 
  coord_polar(theta="y") +
  theme_void()
ggsave(filename = "饼图.pdf", plot = p, device = "pdf", path = "D:/", width = 5, height = 5, dpi = 300, units = "in")
