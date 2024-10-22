# line plot

require(ggplot2)
require(dplyr)

dat <- data.frame(
  time = rep(c(1:4), each =8),
  C = c(rnorm(8, mean = 1.63, sd = 0.05), 
        rnorm(8, mean = 1.81, sd = 0.05), 
        rnorm(8, mean = 1.85, sd = 0.05), 
        rnorm(8, mean = 1.95, sd = 0.05)),
  stringsAsFactors = F
) %>%
  group_by(time) %>%
  mutate(mean = mean(C),
         sd = sd(C)) %>%
  dplyr::select(-C) %>%
  filter(!duplicated(time)) %>%
  mutate(uci = mean + sd,
         lci = mean - sd) %>%
  ungroup()
  
ggplot(data = dat, aes(x = time, y = mean, col = time)) +
  geom_point(shape = 21, size = 2, col = "black") +
  geom_line() +
  geom_errorbar(aes(x = time, ymin = lci, ymax = uci), width = 0.05) +
  xlab("") + 
  ylab("Culture Concentration (nM)") +  
  scale_x_continuous(labels = c("12 h","24 h","48 h","72 h")) +
  scale_y_continuous(limits = c(1.5,2.1), n.breaks = 7, expand = c(0,0)) +
  scale_color_gradient(low = "#b7fce3", high = "#08a269") +
  ggtitle("TMA") +
  guides(size = "none", color = "none") + 
  theme_classic() +
  theme(text = element_text(size = 20), axis.title.y = element_text(vjust = 3),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        plot.margin = ggplot2::margin(50,30,30,30), 
        panel.grid = element_blank(), legend.background = element_blank()
  )

