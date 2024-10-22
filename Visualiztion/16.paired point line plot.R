# paired point line plot

toydat <- data.frame(
  ID = c(1:19),
  `Pre-treatment` = runif(19, min = 1.8, max = 3.6),
  `Post-treatment` = runif(19, min = 0.5, max = 2.8),
 stringsAsFactors = F,
 check.names = F
)

toydat <- toydat %>% 
  reshape2::melt("ID")

paired.t.test <- function(x, y){
  res <- t.test(x, y, alternative = "two.sided", var.equal = TRUE, paired = TRUE)$p.value
  return(res)
  }

# point
ggplot(data = toydat, aes(x = variable, y = value)) +
  geom_line(aes(group = ID), color = "grey50") +
  geom_point(size = 2, color = "#f29b76") + 
  ggsignif::geom_signif(
    test = "t.test",
    test.args = list(paired = TRUE, 
                     alternative = "two.sided",
                     var.equal = F),
    comparisons = list(c("Pre-treatment", "Post-treatment")),
    textsize = 6,
    map_signif_level = TRUE,
    y_position = 3.8
  ) +
  # annotate(geom = "text", x = 1.5, y = 4, label = "***", size = 10) +
  # annotate(geom = "segment", x = 1, xend = 2, y = 3.8, yend = 3.8) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  xlab("") + 
  ylab("Serum testosterone (nmol/L)") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.title = element_blank(),
        plot.margin = ggplot2::margin(30,30,10,30), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(), legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4)
  )

# jitter
ggplot(data = toydat, aes(x = variable, y = value)) +
  geom_line(aes(group = ID), position = position_dodge2(width = 0.2), color = "grey50") +
  geom_point(position = position_dodge2(width = 0.2), size = 2, color = "#f29b76") + 
  ggsignif::geom_signif(
    test = "t.test",
    test.args = list(paired = TRUE, 
                     alternative = "two.sided",
                     var.equal = F),
    comparisons = list(c("Pre-treatment", "Post-treatment")),
    map_signif_level = TRUE,
    textsize = 6,
    y_position = 3.8
  ) +
  scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
  xlab("") + 
  ylab("Serum testosterone (nmol/L)") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.title = element_blank(),
        plot.margin = ggplot2::margin(30,30,10,30), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(), legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4)
  )
