require(plotrix)

ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E',
                 '#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E',
                 '#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E',
                 '#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')

flower_plot <- function(sample, 
                        otu_num, 
                        core_otu, 
                        start, 
                        a, 
                        b, 
                        r, 
                        ellipse_col, 
                        circle_col) {
  
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n  <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1)
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1)
    }
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}