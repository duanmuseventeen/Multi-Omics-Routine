# correlation matrix plot

require(dplyr)
require(GGally)

dat1 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "1") 
dat2 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "2") 
dat3 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "3") 
dat4 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "4") 
dat5 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "5") 

ggpairs(dat1, columns = 2:8, title = "Technical replicates: CLGS")
ggpairs(dat2, columns = 2:8, title = "Technical replicates: IFRV")
ggpairs(dat3, columns = 2:8, title = "Technical replicates: VGTE")
ggpairs(dat4, columns = 2:8, title = "Technical replicates: MVRN")
ggpairs(dat5, columns = 2:8, title = "Technical replicates: NREC")
