# OPLS-DA for each two groups with ropls pkg

# dat with rowname is a df, in which row means sample and column means feature
# group should be factor
opls.each2g <- function(dat, group, output.path = "D:/"){
  suppressMessages(require(ropls))
  suppressMessages(require(ggplot2))
  
  dat.dea <- dat
  group <- as.numeric(group)
  iter <- combn(group %>% unique, 2) %>% as.data.frame
  
  for (i in 1:ncol(iter)) {
    tmp <- dat.dea[group$group %in% iter[[i]],]
    tmp.opls <- opls(
      tmp, 
      group$group[group$group %in% iter[[i]]], 
      permI = 200, 
      scaleC = "pareto",
      crossvalI = 7, 
      predI = 1, 
      orthoI = NA)
    
    if(nrow(tmp.opls@orthoScoreMN) != 0){
      p <- data.frame(
        O1 = tmp.opls@orthoScoreMN[,1],
        PC1 = tmp.opls@scoreMN[,1],
        group = group$group[group$group %in% iter[[i]]] %>% factor,
        label = rownames(tmp),
        row.names = rownames(tmp)) %>% 
        ggplot(aes(PC1, O1, fill = group, color = group)) +
        geom_hline(yintercept = 0, color = "#999999", linetype = 2) + 
        geom_vline(xintercept = 0, color = "#999999", linetype = 2) +
        stat_ellipse(size = 1, geom = "polygon", alpha = 0.1) +
        geom_point(size = 4, alpha = 1) + 
        scale_fill_manual(values = c("#65c294","#b36d41")) +
        scale_color_manual(values = c("#65c294","#b36d41")) +
        guides(fill = "none") +
        labs(x = paste0("t1 (N = ",tmp.opls@descriptionMC[1,1],")\n\n",
                        "R2X = ", tmp.opls@modelDF$`R2X(cum)`[3],"  ",
                        "R2Y = ", tmp.opls@modelDF$`R2Y(cum)`[3],"  ",
                        "Q2Y = ", tmp.opls@modelDF$`Q2(cum)`[3]),
             y = "O1") +
        theme_bw() +
        theme(text = element_text(size = 16))
      ggsave(filename = paste0("OPLS_",paste0(iter[[i]], collapse = "-"),".pdf"), 
             path = output.path, 
             device = "pdf", width = 5, height = 5, dpi = 300)
    }
  }
}
