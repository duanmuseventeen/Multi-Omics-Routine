# https://www.nature.com/articles/s41564-023-01355-5#Sec34
# PMID: 37069399

#install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)

sankey.ed <- readxl::read_excel("4.37069399 Fig5D.xlsx")

sankey_colors<-c("#0072B2", "#999999","#D55E00", "#E69F00","#009E73",  "#56B4E9",  "#F0E442")

ggplot(data = sankey.ed,
       aes(axis1 = `Microbial features`, 
           axis2 = Metabolites, 
           axis3 = `EDI-3 scores`)) +
  scale_x_discrete(limits = c("Microbial features", "Metabolites", "EDI-3 scores")) +
  geom_alluvium(aes(fill = `EDI-3 scores`, alpha =0.8),
                curve_type = "arctangent") +
  geom_stratum(alpha = 0,
               color = "white",
               size=1,
               fill = ) +
  geom_text(stat = "stratum",cex=3, aes(label = after_stat(stratum))) +
  scale_color_manual(values = sankey_colors) +
  scale_fill_manual(values = sankey_colors) +
  guides(fill = "none", alpha = "none") +
  theme_void()+
  theme(legend.position="none",
        axis.text = element_text(size = 16),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid=element_blank())
