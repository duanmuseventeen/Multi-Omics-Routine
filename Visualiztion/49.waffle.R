# 华夫饼图----
# https://github.com/hrbrmstr/waffle
# devtools::install_github("hrbrmstr/waffle")
library(waffle)
library(magrittr)
library(hrbrthemes)

dat <- readxl::read_excel("metabolite classified by HMDB.xlsx")

head(dat)
# A tibble: 6 × 8
# Query                      `HMBD Class`
# <chr>                       <chr>       
#   1 1-Methylnicotinamide   Pyridines a…
# 2 3-Hydroxybutyric acid    Beta hydrox…

dat %>%
  mutate(vals = 1) %>% 
  count(`HMBD Class`, wt = vals) %>%
  ggplot(aes(fill = `HMBD Class`, values = n)) +
  geom_waffle(
    n_rows = 13,
    size = 0.33, 
    colour = "white",
    flip = TRUE) +
  coord_equal() +
  theme_ipsum_rc(grid="") +
  theme_enhance_waffle()
