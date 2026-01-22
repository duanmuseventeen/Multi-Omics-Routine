# https://www.zhihu.com/question/305623922

# gradient ------------------------------------------------------------------------
colorRampPalette(c("white","#007947"))(200)
                
colorRampPalette(c("white","red3"))(200)        

# (PMID: 40345201)
c("black","red","orange")
c("black","blue","green")

# warm (PMID: 38458196)
colorRampPalette(
    c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
      '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
      '#6F003D','#56033F'))(1000)

# (PMID: 38959864)
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

# (PMID: 40670783)
c(low = "#8b81bd", mid = "white", high = "#822425")

# (PMID: 40301344)
 scale_fill_distiller(palette = "RdBu",
                      limits = c(-0.8, 0.8),
                      breaks = seq(-0.8, 0.8, 0.4),
                      direction = 1) 

# for single cell -----------------------------------------------------------------
# https://github.com/friedpine/scRNASeq_ESCC/blob/main/3c.R
mycol <-c("#DEEAB2","#64B473","#2D553C","#A1D8C9","#487C76","#7AAF93","#D0E4C6",
          "#F3746C","#BB4B95","#F8BFAF","#F7DDD3","#66CDF6","#598198","#D5E7F7",
          "#69B3CE","#D6D5EB","#7B8CBC","#7674AE", "#E3CEE4",
          "#AFB2B6","#C9BDB2","#F5C785","#ECAECF","#E9A943","#CAA57D","#A79388",
          "#EACF68","#F6F3A7","#C45337","#86382B","#EADCE4",
          "#EE5276","#9E6CAB","#74507B")

# https://mp.weixin.qq.com/s?__biz=MzkwNTY1ODMwMw==&mid=2247487820&idx=1&sn=97bff731a76911156492ca40521d0ef9&chksm=c0f53b9cf782b28ac865935db146048ebe0146035f97f388672f5597a0f7ef6021ea9681e5c3&scene=21#wechat_redirect
c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442",
  "#0072B2","#D55E00","#CC79A7","#CC6666",
  "#9999CC","#66CC99","#3C5488B2","#00A087B2", 
  "#F39B7FB2","#91D1C2B2", 
  "#8491B4B2", "#DC0000B2", 
  "#7E6148B2","yellow", 
  "darkolivegreen1", "lightskyblue")

c("#e8d5d5","#fad8b7","#f5b0af","#b3ccb0","#d2e7a7","#afe5f6")
c("#c79494","#f39d49","#e73a36","#3d7d35","#8ec321","#30bbea")

c("#F37252","#F7935A","#FCB461","#FED297","#FFF0DC","#E3F5F6","#ABE0E4","#7FC6D3","#779FC6","#6F78B9")


new_colors_list = {
 'IPHC/IBC' : '#E76254',
 'Fibroblasts': '#EF8A47',
 'DC/PC': '#f4a494',
 'Basal cell': '#FFE6B7',
 'Intermediate cell': '#AADCE0',
 'Marginal cell': '#528FAD',
 'Hair cell': '#a4549c',
 'Spiral ganglion neuron': '#1E466E',
 'Reissner’s membrane': '#C7B8BD',
 'Schwann cell': '#8C4834',
 'Capillary endothelial cell': '#C17E65',
 'Pericytes': '#645cac',
 'T/NK': '#EFD061',
 'B cell': '#547857',
 'Mac/Mono': '#c49c94',
 'Red blood cell': '#f7b6d2',
 'Neutrophil': '#dbdb8d'
}





