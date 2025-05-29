# https://wilkelab.org/ggridges/articles/introduction.html

require(ggridges)

iris |>
  ggplot(aes(x = Sepal.Length, y = Species, fill = Species)) + 
  geom_density_ridges(scale = 1) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")))

iris |>
  ggplot(aes(x = Sepal.Length, y = Species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1) +
  scale_fill_viridis_c(name = "Temp. [F]", option = "C")

