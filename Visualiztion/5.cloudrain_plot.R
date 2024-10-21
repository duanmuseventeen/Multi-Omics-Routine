require(ggExtra)
require(ggplot2)
require(ggpubr)
require(dplyr)
require(stringr)
require(gghalves)
require(gglayer)
require(PupillometryR)

# method 1
p1 <- ggplot(iris, aes(Species, Petal.Length, color = Species, fill = Species)) +
  geom_flat_violin(position = position_nudge(x = .3), color = NA, alpha = 0.2) +
  geom_jitter(aes(color = Species), width = .15, alpha = 0.2) +
  geom_boxplot(width = .1, position = position_nudge(x = .22), alpha = 0.2) +
  labs(y = "Shannon Index", x = "") +
  guides(color = "none", fill = "none") +
  coord_flip() +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    plot.margin = ggplot2::margin(40, 40, 40, 40),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(plot = p1, filename = "fig1b.pdf", path = "D:/", device = "pdf", width = 6, height = 6, dpi = 300)

# method 2
require(raincloudplots)
require(ggsci)
require(ggtext)
df_1x1 <- data_1x1(
  array_1 = iris$Sepal.Length[1:50],
  array_2 = iris$Sepal.Length[51:100],
  jit_distance = .09,
  jit_seed = 321
)
mycolor <- pal_jco()(2)
raincloud_1_h <- raincloud_1x1(
  data = df_1x1,
  colors = mycolor,
  fills = mycolor,
  size = 1,
  alpha = .8,
  ort = "h"
)
raincloud_2 <- raincloud_1x1_repmes(
  data = df_1x1,
  colors = mycolor,
  fills = mycolor,
  line_color = "gray",
  line_alpha = .6,
  size = 1,
  alpha = .8,
  align_clouds = FALSE
) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Pre", "Post"), limits = c(0, 3)) +
  xlab("Time") +
  ylab("Score") +
  labs(
    title = "Example of <span style='color:#D20F26'>raincloudplots::raincloud_1x1 function</span>",
    subtitle = "processed charts with <span style='color:#1A73E8'>raincloud_1x1()</span>",
    caption = "Visualization by <span style='color:#0057FF'>DataCharm</span>"
  ) +
  hrbrthemes::theme_ipsum(base_family = "Roboto Condensed") +
  theme(
    plot.title = element_markdown(
      hjust = 0.5, vjust = .5, color = "black",
      size = 20, margin = margin(t = 1, b = 12)
    ),
    plot.subtitle = element_markdown(hjust = 0, vjust = .5, size = 15),
    plot.caption = element_markdown(face = "bold", size = 12),
    legend.position = "none"
  )
