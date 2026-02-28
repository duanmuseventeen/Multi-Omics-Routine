# chord plot
rm(list = ls());gc()

# Usually, this plot are performed from a sequencing data processed by limma, DESeq2
# or edgeR. Here, I just use GOChord to perform a chord plot to display the association
# between genes and pathways. Of course, you can use this way to display any
# relationship between category and units, for example metabolites.

# load pkgs---------------------------------------------------------------------
require(ggtext)
require(dplyr)
require(clusterProfiler)
require(enrichplot)
require(GOplot)
require(stringr)

# graph-------------------------------------------------------------------------
# I create the pathway and genelist files according to tutorial. For our purpose
# some information are not necessary, they are left empty
pathway <- readxl::read_excel("pathway.xlsx")
genelist <- readxl::read_excel("genelist.xlsx")

circ <- circle_dat(pathway, genelist)

chord <- chord_dat(data = circ, genes = genelist %>% mutate(ID = toupper(ID)) %>% as.data.frame)

pdf("chord.pdf", width = 9, height = 10)
mygochord(chord, 
        space = 0.02, 
        gene.order = 'logFC', 
        gene.space = 0.25, 
        gene.size = 5,
        ribbon.col = c("#F37252", "#F7935A", "#FCB461", "#FED297", "#FFF0DC",
                       "#E3F5F6", "#ABE0E4", "#7FC6D3", "#779FC6"))
dev.off()

# for ggplot2 update-------------------------------------------------------------
check_chord <- function(mat, limit){
  
  if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
  
  tmp <- mat[(rowSums(mat) >= limit[1]),]
  mat <- tmp[,(colSums(tmp) >= limit[2])]
  
  mat <- check_chord(mat, limit)
  return(mat)
}

bezier <- function(data, process.col){
  x <- c()
  y <- c()
  Id <- c()
  sequ <- seq(0, 1, by = 0.01)
  N <- dim(data)[1]
  sN <- seq(1, N, by = 2)
  if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
  for (n in sN){
    xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
    for (t in sequ){
      xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
      xval <- c(xval, xva)
      xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
      xval2 <- c(xval2, xva2)
      yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]  
      yval <- c(yval, yva)
      yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
      yval2 <- c(yval2, yva2)			
    }
    x <- c(x, xval, rev(xval2))
    y <- c(y, yval, rev(yval2))
    Id <- c(Id, rep(n, 2 * length(sequ)))
  }
  df <- data.frame(lx = x, ly = y, ID = Id)
  return(df)
}

mygochord <- function (data, title, space, gene.order, gene.size, gene.space, 
          nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, 
          process.label, limit) 
{
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  if (missing(title)) 
    title <- ""
  if (missing(space)) 
    space = 0
  if (missing(gene.order)) 
    gene.order <- "none"
  if (missing(gene.size)) 
    gene.size <- 3
  if (missing(gene.space)) 
    gene.space <- 0.2
  if (missing(lfc.col)) 
    lfc.col <- c("brown1", "azure", "cornflowerblue")
  if (missing(lfc.min)) 
    lfc.min <- -3
  if (missing(lfc.max)) 
    lfc.max <- 3
  if (missing(border.size)) 
    border.size <- 0.5
  if (missing(process.label)) 
    process.label <- 11
  if (missing(limit)) 
    limit <- c(0, 0)
  if (gene.order == "logFC") 
    data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == "alphabetical") 
    data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
    if (nlfc == 1) {
      cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x, 
                                                            rownames(data)), Ncol])
    }
    else {
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, 
                                                      (Ncol - nlfc + 1)])
    }
  }
  else {
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col)) 
    colRib <- grDevices::rainbow(dim(cdata)[2])
  else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 
                                                      202 * nrib[b]))
  r1 <- 1
  r2 <- r1 + 0.1
  xmax <- c()
  x <- 0
  for (r in 1:length(nrib)) {
    perc <- nrib[r]/sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) 
      x <- c(x, x[r] + pi * perc)
  }
  xp <- c()
  yp <- c()
  l <- 50
  for (s in 1:Ncol) {
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + 
                                                         xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), 
            r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + 
                                                         xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), 
            r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), 
                                                    each = 4 + 2 * l))
  xp <- c()
  yp <- c()
  logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi/Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi/Nrow + space, length = Nrow)
  for (s in 1:Nrow) {
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1) {
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * 
                sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), 
              r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), 
              r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * 
                cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), 
              r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), 
              r2 * cos(x2[s]))
    }
    else {
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc) {
        logs <- c(logs, data[s, (dim(data)[2] + 1 - 
                                   t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * 
                  sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), 
                tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * 
                  cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), 
                tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }
    }
  }
  if (lfc[1] != 0) {
    if (nlfc == 1) {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                      each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 
                                                                                       2 * l))
    }
    else {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc * 
                                                             Nrow)), each = 4 + 2 * l), logFC = rep(logs, 
                                                                                                    each = 4 + 2 * l))
    }
  }
  else {
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                    each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2))
  angle <- c()
  for (o in aseq) if ((o + 270) <= 360) 
    angle <- c(angle, o + 270)
  else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + 
                                                         xmax2/2), ygen = (r1 + gene.space) * cos(x2 + xmax2/2), 
                        labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax/2), 
                        ypro = (r1 + 0.15) * cos(x + xmax/2), labels = colnames(cdata), 
                        stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c()
  y.end <- c()
  processID <- c()
  for (gs in 1:length(x2)) {
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 
                 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)) {
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), 
  ]
  x.start <- c()
  y.start <- c()
  for (rs in 1:length(x)) {
    val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 
                 1)
    for (v in 1:(length(val) - 1)) {
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if (length(df_genes$logFC) != 0) {
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x > 
                                                       lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, 
                                            lfc.min, x))
    df_genes$logFC <- logFC
  }
  g <- ggplot() + geom_polygon(data = df_process, aes(x, y, 
                                                      group = id), fill = "gray70", inherit.aes = F, color = "black") + 
    geom_polygon(data = df_process, aes(x, y, group = id), 
                 fill = cols, inherit.aes = F, alpha = 0.6, color = "black") + 
    geom_point(aes(x = xpro, y = ypro, size = factor(labels, 
                                                     levels = labels), shape = NA), data = df_texp) + 
    guides(size = guide_legend("GO Terms", ncol = 4, byrow = T, 
                               override.aes = list(shape = 22, fill = unique(cols), 
                                                   size = 8))) + theme(legend.text = element_text(size = process.label)) + 
    geom_text(aes(xgen, ygen, label = labels, angle = angle), 
              data = df_texg, size = gene.size) + geom_polygon(aes(x = lx, 
                                                                   y = ly, group = ID), data = df_path, fill = colRibb, 
                                                               color = "black", size = border.size, inherit.aes = F) + 
    labs(title = title) + theme_void()
  if (nlfc >= 1) {
    g + geom_polygon(data = df_genes, aes(x, y, group = id, 
                                          fill = logFC), inherit.aes = F, color = "black") + 
      scale_fill_gradient2("logFC", space = "Lab", low = lfc.col[3], 
                           mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", 
                                                                                       title.hjust = 0.5), breaks = c(min(df_genes$logFC), 
                                                                                                                      max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), 
                                                                                                                                                       round(max(df_genes$logFC)))) + theme(legend.position = "bottom", 
                                                                                                                                                                                            legend.background = element_rect(fill = "transparent"), 
                                                                                                                                                                                            legend.box = "horizontal", legend.direction = "horizontal")
  }
  else {
    g + geom_polygon(data = df_genes, aes(x, y, group = id), 
                     fill = "gray50", inherit.aes = F, color = "black") + 
      theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent"), 
            legend.box = "horizontal", legend.direction = "horizontal")
  }
}

