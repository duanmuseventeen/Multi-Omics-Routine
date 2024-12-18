require(ggtext)
require(ggplot2)
require(dplyr)

p <- readxl::read_excel("***.xlsx") %>% 
  mutate(sig = case_when(FC >= log2(1.5) & `p-value` < 0.05 ~ 2,
                         FC <=-log2(1.5) & `p-value` < 0.05 ~ 1,
                         TRUE ~ 0),
         sig = factor(sig, labels = c("Non-sig","Down","Up"))) %>% 
  ggplot(aes(x = FC, y = -log10(`p-value`), col = sig, fill = sig)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "gray50", linetype = 2) +
  geom_hline(yintercept = -log10(0.05), col = "gray50", linetype = 2) +
  geom_point(size = 2, shape = 21, alpha = 0.5) +
  scale_color_manual(values = c("gray70","blue3","red3"), name = "") +
  scale_fill_manual(values = c("gray70","blue3","red3"), name = "") +
  labs(x = "log<sub>2</sub>(Fold Change)<br>Model vs Control", 
       y = "-log<sub>10</sub>(<i>p</i>value)") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_text(hjust = 0.5, vjust = 5),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.background = element_blank(),
        plot.margin = ggplot2::margin(50,30,30,30), 
        legend.position = "top",
        legend.background = element_blank())

ggsave(plot = p, filename = "***.pdf", path = "", device = "pdf", width = 5, height = 5, dpi = 300)