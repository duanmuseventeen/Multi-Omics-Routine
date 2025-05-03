require(ggplot2)
require(ggnewscale)
require(dplyr)

dat <- readxl::read_excel("7.slingshot.xlsx")

p <- ggplot() +
  geom_point(data = dat %>% filter(complete.cases(slingPseudotime_1)),
             aes(x = umap_1, y = umap_2, col = slingPseudotime_1, alpha = 0.1)) +
  scale_color_viridis(name="Slingshot:\nPseudotime1", option = "D") +
  ggnewscale::new_scale_color() +
  geom_point(data = dat %>% filter(complete.cases(slingPseudotime_2)),
             aes(x = umap_1, y = umap_2, col = slingPseudotime_2, alpha = 0.1)) +
  geom_path(data = dat.fit1$s[dat.fit1$ord,],
            aes(x = umap_1, y = umap_2), col = "red3") +
  geom_path(data = dat.fit2$s[dat.fit2$ord,],
            aes(x = umap_1, y = umap_2), col = "blue3") +
  scale_color_viridis(name="Slingshot:\nPseudotime2", option = "A") +
  ggtitle("Slingshot CD8T") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = ggplot2::margin(15,15,15,15),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(filename = "CD8T cell(Slingshot).pdf", plot = p, device = "pdf",
       width = 7, height = 5, dpi = 300, units = "in")
