dat <- data.frame(
  value = runif(40, min = 0, max = 200),
  group = c(rep("sham",10), rep("model", 10), rep("low", 10), rep("high",10)),
  condition = c(rep(paste0("sample_",c(1:10)), 4)),
  log2FC = runif(40, min = -2, max = 2),
  p = runif(40, min = 0, max = 0.25)
  ) %>% 
  mutate(size = -log10(p),
         sig = ifelse(p > .05, "", 
                      ifelse(p > 0.01, "*", 
                             ifelse(p > 0.001, "**", "***"))))

ggplot(dat, aes(group, condition, color = log2FC, size = size, label = sig)) +
  geom_point() +
  geom_text(color = "#222222", size = 3) +
  labs(x = "", y = "") +
  scale_size_continuous(range = c(1,6)) +
  scale_color_gradientn(
    breaks = seq(-2,2,0.5),
    limits = c(-2,2),
    name = "",
    colors = colorRampPalette(c("blue4", "white", "red3"))(9)) +
  theme_minimal() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p, filename = "bubble.pdf", path = "F:/", 
       device = "pdf", width = 10, height = 10, dpi = 300)