require(ggplot2)
require(dplyr)

# Ref: 
# https://blog.csdn.net/weixin_43151909/article/details/106450585
# example 1---------------------------------------------------------------------
ggplot(data = mtcars, aes(x = mpg, y = wt)) + 
  geom_point() + 
  geom_polygon(data = mtcars[chull(mtcars$mpg,mtcars$wt),], aes(x = mpg, y = wt),
               alpha = 0.2)

# example 2---------------------------------------------------------------------
mychull <- function(dat) {dat[chull(dat$mpg, dat$wt),]}

ggplot(data = mtcars, aes(x = mpg, y = wt, color= factor(cyl), fill = factor(cyl))) + 
  geom_point() + 
  geom_polygon(data = plyr::ddply(mtcars, ~ cyl, mychull),
               aes(x = mpg, y = wt, color= factor(cyl)),
               alpha = 0.2) + 
  labs(color = "", fill = "") +
  theme_bw() +
  theme(plot.margin = margin(c(30,30,30,30)), panel.grid = element_blank())

ggplot(data = mtcars, aes(x = mpg, y = wt, color= factor(cyl), fill = factor(cyl))) + 
  geom_point() + 
  geom_polygon(data = mtcars %>%
                 group_by(cyl) %>%
                 slice(chull(mpg, wt)),
               aes(x = mpg, y = wt, color= factor(cyl)),
               alpha = 0.2) + 
  labs(color = "", fill = "") +
  theme_bw() +
  theme(plot.margin = margin(c(30,30,30,30)), panel.grid = element_blank())










