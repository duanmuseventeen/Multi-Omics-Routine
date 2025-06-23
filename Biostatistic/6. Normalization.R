# https://www.jianshu.com/p/de25afe02a33
# https://www.cnblogs.com/yanzhi123/p/11712926.html
# https://www.metaboprofile.com/school/1935/
# pmid: 24766612
# pmid: 24942700
# pmid: 27694351

# normalization-----------------------------------------------------------------
for (i in 1:ncol(data)) {
  tmp <- data[[i]]
  
  # mean normalization 
  data[[i]] <- tmp / mean(tmp)
  # median normalization 
  data[[i]] <- tmp / median(tmp)
  # max normalization 
  data[[i]] <- tmp / max(tmp)
  # global normalization (归一化)
  data[[i]] <- tmp / sum(tmp)
  # max-min normalization 
  data[[i]] <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
  # Z-score normaliztion
  data[[i]] <- (tmp - mean(tmp)) / sd(tmp)
  # IQR normaliztion
  data[[i]] <- (tmp - median(tmp)) / iqr(tmp)
}

# VSN
limma::normalizeVSN()

# scale-------------------------------------------------------------------------
# Center scaling
x - mean(x)

# unit variance scaling
(x - mean(x)) / sd(x)
scale(x)

# Pareto scaling
(x - mean(x)) / sqrt(sd(x))

# scale to 0 ~ 1
(x - min(x)) / (max(x) - min(x))

# scale to -1 ~ 1
(x - min(x)) / (max(x) - min(x)) * 2 - 1











