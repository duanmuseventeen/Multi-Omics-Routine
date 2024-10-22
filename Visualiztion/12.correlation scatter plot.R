# correlation scatter plot

require(dplyr)
require(tibble)
require(ggplot2)
require(ggtext)  
require(ggExtra)

dat <- tibble::tibble(
  x = runif(16, min = 17, max = 78),
  y = runif(16, min = 0.4, max = 1.15)
)

# eq <- substitute(expr = paste("Relative ", italic("P"), ".copri expression level")) #must use paste()

p <- ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, col = "black", fill = "#94c4fd") +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, col = "black") + 
  xlab("Plasma TMAO level (μM)") + 
  # method 1
  # ylab(eq) +
  # method 2 # require(ggtext)
  ylab("Relative <i>P</i>. copri expression level") +  
  geom_richtext(
    data = data.frame(
      x = 30, 
      y = 0.2, 
      label = "<b>Spearman r = 0.7529<br/>n=16 <i>p</i>=0.001</b>"),
    aes(x, y, label = label), 
    size = 6, 
    hjust = 0,
    fill = NA, 
    label.color=NA, # remove bubble
    show.legend = F) +
  scale_x_continuous(limits = c(0,80), breaks = c(0,20,40,60,80), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.5), breaks = c(0,0.5,1,1.5), expand = c(0,0)) +
  ggtitle("FMT") +
  guides(size = "none", color = "none") + 
  # coord_cartesian(clip = "off") + 
  theme_classic() +
  theme(text = element_text(size = 20), 
        legend.title = element_blank(),
        axis.title.y = element_markdown(),
        plot.title = element_text(hjust = 0.5, vjust = 5),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.background = element_blank(),
        plot.margin = ggplot2::margin(50,30,30,30), 
        panel.grid = element_blank(), 
        legend.background = element_blank()
  )
  
ggMarginal(p, 
           type = "densigram"
           # groupColour = TRUE,
           # groupFill = TRUE
           )