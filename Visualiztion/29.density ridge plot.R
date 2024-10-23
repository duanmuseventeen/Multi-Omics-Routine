# ridge plot

require(ggridges)

iris |>
  ggplot(aes(x = Sepal.Length, y = Species, fill = Species)) + 
  geom_density_ridges(scale = 1)
