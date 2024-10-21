library(venn)

sets <- as.data.frame(matrix(sample(0:1, 150, replace = TRUE), ncol = 5))
colnames(sets) <- c("Glucose 6G", "Glucose 15G", "Leucine", "Fatty acid", "KCl")
  
warmcol <- c("#FF5733", #(a deep orange-red)
             "#B22222", #"firebrick",
             #"#FFA07A", #(a light salmon pink)
             #"#FFDAB9", #(a peachy beige)
             "#FFC0CB", #"pink"
             "#FFC300", #(a bright yellow-orange)
             "#FF8C00") #(a bright orange)

coldcol <- c(#"turquoise2",
  #"#66CDAA", #"aquamarine3"
  #"#7FFFD4",#"aquamarine",
  #"#00BFFF",#"medium blue",
  #"#4169E1",#"royal blue",
  "#4682B4",#, "steelblue",
  "#6A5ACD",# "slateblue,
  #"#7FFFD4",#"light turquoise",
  "#B0E0E6",# "powderblue",
  #"#00FFFF",#"bright cyan", 
  #"#E0FFFF", #"pale blue-green"
  #"#6a51a3", #medianpurple4
  "#5F9EA0", #"cadetblue"
  "skyblue2"
)

venn(sets, 
     box = F, 
     ilabels = "counts",
     #ilab=TRUE, 
     snames = c("Glucose 6G", "Glucose 15G", "Leucine", "Fatty acid", "KCl"),
     zcolor = warmcol, #"style", #mycol,
     borders = T,
     ilcs = 2, # control size of intersection numbers
     sncs = 2, # control size of group names
     ggplot = T,
     size=1, 
     color=warmcol #"tomato",
)

ggsave(filename = "venn.pdf", width=10, height=10, units="in")