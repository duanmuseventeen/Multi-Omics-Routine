# Bubble plot

require(ggplot2)
require(dplyr)

dat <- readxl::read_excel("8-38347143-Source Data Extended Fig8.xlsx",
                          sheet = "ED_8c") %>%
  arrange(`Importance Score`) %>% 
  mutate(`Feature type` = rev(c("Clinical stage", rep("Metabolite", 4), "IHC subtype",
                            rep("RNA", 2), rep("Pathology", 3), rep("RNA", 2),
                            "Pathology", "Metabolite", "RNA", "Pathology")),
         Feature = factor(Feature, levels = Feature)) %>% 
  mutate(col = case_when(`Feature type` == "Clinical stage" ~ "yellow3",
                         `Feature type` == "Metabolite" ~ "red3",
                         `Feature type` == "IHC subtype" ~ "blue3",
                         `Feature type` == "RNA" ~ "green3",
                         `Feature type` == "Pathology" ~ "purple3"))
ggplot(dat, aes(x = `Importance Score`, y = Feature, 
           color = `Feature type`,
           fill = `Importance Score`)) +
  geom_point(size = 5, 
             pch = 21, 
             col = "gray90") +
  scale_fill_gradient(low = "white", high = "#6683a7") +
  labs(x = "Feature Importance Score", y = "") +
  scale_y_discrete() +
  scale_x_continuous(limits = c(0, 0.18), breaks = c(0.05,0.1,0.15), expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.frame = element_rect(color = "black"),
        axis.text.y = element_markdown(color = dat$col),
        plot.margin = ggplot2::margin(30,30,30,30))
