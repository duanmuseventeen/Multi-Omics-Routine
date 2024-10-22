# boxplot

require(ggplot2)
require(dplyr)

dat <- readxl::read_excel("23.38347143-Source Data Extended Fig7.xlsx", sheet = "ED_7d,e")

dat %>% 
  mutate(x = factor(Clinical_Subtype, 
                    levels = c("HR+HER2-","HR+HER2+","HR-HER2+","TNBC","PT"),
                    labels = c("HR+HER2-","HR+HER2+","HR-HER2+","TNBC","Paratumour"))) %>% 
  ggplot(aes(x, `Number of fusion transcripts`, colour = x)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) + 
  scale_color_manual(values = c("#1984be", "#7eaa42", "#dba046","#c9543c","#87add0")) +
  labs(x = "", y = "Fusion transcripts per sample") +
  guides(col = "") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position = c(0.9,0.6))

dat %>% 
  mutate(x = factor(PAM50_classifier, 
                    levels = c("LumA","LumB","Her2","Basal","Normal","PT"),
                    labels = c("Luminal A", "Luminal B", "HER2-enriched",
                               "Basal-like", "Normal-like", "Paratumour"))) %>% 
  ggplot(aes(x, `Number of fusion transcripts`, colour = x)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) + 
  scale_color_manual(values = c("#1976bc", "#64ccdc", "#655ead","#dd223a","#c5c6c8","#87aecd")) +
  labs(x = "", y = "Fusion transcripts per sample") +
  guides(col = "") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position = c(0.9,0.6))
