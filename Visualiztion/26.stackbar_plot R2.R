# barplot

require(ggplot2)
require(tidyr)
require(dplyr)

dat <- readxl::read_excel("25.38347143-Source Data Extended Fig6.xlsx", sheet = "ED_6a")

dat1 <- dat %>% 
  filter(`...1` != "TME phenotype") %>% 
  pivot_longer(cols = ACKR:ZYWC, names_to = "Sample", values_to = "value") %>% 
  transmute(`Cell type` = `...1`, 
            Sample, value) %>% 
  mutate(value = as.numeric(value)) 

dat2 <- dat %>%
  filter(`...1` == "TME phenotype") %>% 
  tibble::column_to_rownames("...1") %>% 
  t %>% as.data.frame %>% 
  tibble::rownames_to_column("Sample")

dat.merge <- dat1 %>% 
  left_join(dat2, by = "Sample") %>% 
  mutate(`TME phenotype` = factor(`TME phenotype`, levels = c("Cold", "Moderate", "Hot"))) %>% 
  arrange(`TME phenotype`) %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample)))

ggplot(dat.merge, aes(Sample, value, fill = `Cell type`)) +
  geom_col(alpha = 0.8) +
  geom_segment(data = dat.merge, 
               aes(x = Sample, xend = Sample,
                   y = 1.05, yend = 1.1, 
                   colour = `TME phenotype`)) +
  labs(x = "", y = "Proportion") +
  guides(fill = guide_legend(ncol = 1)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(text = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 16, angle = 90, vjust = 4), 
        axis.text.y = element_text(size = 12), 
        axis.ticks.y = element_line(linewidth = 1, color = "black"),
        plot.margin = ggplot2::margin(20,20,20,50),
        legend.position = "right", legend.title = element_text(face = "bold"))
  
