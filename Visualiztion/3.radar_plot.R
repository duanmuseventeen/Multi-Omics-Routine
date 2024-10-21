# devtools::install_github("ricardo-bion/ggradar")
library("ggradar")
library(tidyverse)

load("data for radar.Rdata")

ggradar(
  df, 
  values.radar = c("-400%", "0%", "150%"),
  grid.min = 0, grid.mid = 400, grid.max = 550,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
  fill = TRUE,
  fill.alpha = 0.2,
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey"
)

ggsave(filename = paste0("RADAR",".pdf"), 
       path = "D:/", device = "pdf",
       width = 15, height = 8, dpi = 300)
