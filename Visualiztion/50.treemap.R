dat <- readxl::read_excel("metabolite classified by HMDB.xlsx")

# treemapify----
# install.packages("treemapify")
library(treemapify)

dat %>%
  mutate(vals = 1) %>% 
  count(`HMBD Class`, value = vals) %>%
  group_by(`HMBD Class`) %>% 
  mutate(prop = n / nrow(dat)) %>% 
  mutate(prop = format(prop * 100, digits = 2) %>% paste0(.,"%")) %>% 
  ggplot(aes(area = value, fill = `HMBD Class`, 
             label = prop, subgroup = `HMBD Class`)) +
  geom_treemap() +
  geom_treemap_subgroup_border(colour = "white", size = 5) +
  geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.5, colour = "black") +
  geom_treemap_text(colour = "white", place = "topleft", reflow = TRUE) +
  # scale_fill_brewer(palette = "Set2") +
  labs(title = "")
# treemap----
# install.packages("treemap")
library(treemap)

dat %>%
  mutate(vals = 1) %>% 
  count(`HMBD Class`, value = vals) %>%
  group_by(`HMBD Class`) %>% 
  mutate(prop = n / nrow(dat)) %>% 
  mutate(prop = format(prop * 100, digits = 2) %>% paste0(.,"%")) %>% 
  treemap(
    index = "HMBD Class", 
    vSize = "value",                   
    vColor = "HMBD Class",               
    type = "index",                    
    title = "Simple Treemap",
    palette = "Pastel1")
