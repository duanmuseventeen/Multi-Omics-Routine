# 代谢物同异名的问题是代谢组学处理中的一大痛点，该脚本提供一个可行的解决办法，但无法完全解决所有可能出现的问题，更多特殊情况仍需人工校对

library(webchem)
library(dplyr)

dat   <- readxl::read_excel("1.metabolite2align.xlsx")
cid   <- dat$Metabolite %>% unique %>% get_cid
dat   <- dat %>% 
  left_join(cid %>% mutate(Metabolite = query), by = "Metabolite")

export::table2excel(dat, "1.metabolite2align(add cid).xlsx")
# 后续再Excel中根据cid进行代谢物对齐即可
