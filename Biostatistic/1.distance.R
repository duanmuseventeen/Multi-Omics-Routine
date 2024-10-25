# distance

toydat <- data.frame(
  a = runif(5, 0,   1),
  b = runif(5, 0.3, 1.5),
  c = runif(5, 0.5, 1.3)
)

toydat
# 1. Euclidean distance---------------------------------------------------------
toydat |> 
  apply(1, function(x){
    apply(toydat, 1, function(y){
      sum((x - y)^2)^0.5
    })
  })

dist(toydat, method = "euclidean")
# 2. Maximum distance---------------------------------------------------------
toydat |> 
  apply(1, function(x){
    apply(toydat, 1, function(y){
      max(abs(x - y))
    })
  })

dist(toydat, method = "maximum")
# 3. Manhattan distance---------------------------------------------------------
toydat |> 
  apply(1, function(x){
    apply(toydat, 1, function(y){
      sum(abs(x - y))
    })
  })

dist(toydat, method = "manhattan")
# 4. Canberra distance---------------------------------------------------------
toydat |> 
  apply(1, function(x){
    apply(toydat, 1, function(y){
      sum(abs(x - y) / (abs(x) + abs(y)))
    })
  })

dist(toydat, method = "canberra")
# 5. Minkowski distance---------------------------------------------------------
p <- 3
toydat |> 
  apply(1, function(x){
    apply(toydat, 1, function(y){
      sum(abs(x - y) ^ p)^(1/p)
    })
  })

dist(toydat, method = "minkowski", p = 3)






