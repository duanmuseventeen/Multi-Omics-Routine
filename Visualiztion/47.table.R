# Method1-----------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/134034762

mytable <- mtcars %>% 
  group_by(cyl) %>% 
  summarise(sum(mpg))

ggplot(data = mtcars, aes(x = factor(cyl), y = mpg, fill = factor(cyl))) +
  geom_col() +
  annotation_custom(grob = gridExtra::tableGrob(mytable, rows = NULL),
                    xmin = 7, xmax = 9,
                    ymin = 250, ymax = 300) +
  theme_bw() +
  theme(text = element_text(size = 14))
# Method2-----------------------------------------------------------------------
# https://github.com/aphalo/ggpp
library(ggpp)

mtcars %>%
  group_by(cyl) %>%
  summarize(wt = mean(wt), mpg = mean(mpg)) %>%
  ungroup() %>%
  mutate(wt = sprintf("%.2f", wt),
         mpg = sprintf("%.1f", mpg)) -> tb

df <- tibble(x = 5.45, y = 34, tb = list(tb))

ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point() +
  geom_table(data = df, aes(x = x, y = y, label = tb))
