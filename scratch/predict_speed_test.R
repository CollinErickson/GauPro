# Create model
n <- 40
d <- 4
x <- matrix(runif(n*d), ncol=d)
f1 <- function(a) {sin(2*pi*a[1]) + sin(6*pi*a[2]) + abs(sin(2*pi*a[3])) + a[4]}
y <- apply(x,1,f1)# + rnorm(n,0,.01)
gp <- GauPro(x,y, verbose=2);gp$theta


N <- 1e4
# n <- c(1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,300,1e3)
n <- c(30,35,40,45,50,55,1e3)
Nn.floor <- floor(N/n)
times <- rep(0, length(n))


# Create prediction data
XX <- matrix(runif(N*d), ncol=d)
for (i in seq_along(n)) {
  st <- system.time({
  ni <- n[i]
  Nni <- ceiling(N/ni)-1
  #sapply(1:Nni, function(j) {print(c(ni,(j*ni+1), (min((j+1)*ni,N))))})
  #tmp <- sapply(1:Nni, function(j) {(gp$predict(XX[(j*ni+1):(min((j+1)*ni,N)), ]))})
  predout <- numeric(N)
  for (j in 0:Nni) {
    predout[(j*ni+1):(min((j+1)*ni,N))] <- gp$pred(XX[(j*ni+1):(min((j+1)*ni,N)), ], split_speed = F)
  }
  #sapply(1:Nni, function(j) {predout[(j*ni+1):(min((j+1)*ni,N))] <<- (gp$predict(XX[(j*ni+1):(min((j+1)*ni,N)), ]))})
  })
  cat(paste(i, ni, st[1], '\n'))
}


# For 2D data
# i, ni, st[1]
# [1] "1 1 5.96"
# [1] "2 3 2.41"
# [1] "3 10 1.03"
# [1] "4 30 0.609999999999999"
# [1] "5 100 0.859999999999999"
# [1] "6 300 1.41"
# [1] "7 1000 4.26"
# [1] "8 3000 12.13"
# [1] "9 10000 32.23"



# [1] "1 1 5.67"
# [1] "2 2 2.53999999999999"
# [1] "3 3 1.7"
# [1] "4 5 1.13000000000001"
# [1] "5 10 0.75"
# [1] "6 20 0.509999999999991"
# [1] "7 30 0.560000000000002"
# [1] "8 40 0.509999999999991"
# [1] "9 50 0.569999999999993"
# [1] "10 60 0.489999999999995"
# [1] "11 70 0.52000000000001"
# [1] "12 80 0.469999999999999"
# [1] "13 90 0.609999999999999"
# [1] "14 100 0.61999999999999"
# [1] "15 300 1.08"
# [1] "16 1000 3.33"


# [1] "1 10 6.87"
# [1] "2 20 5"
# [1] "3 30 4.60000000000001"
# [1] "4 40 4.39"
# [1] "5 50 4.47999999999999"
# [1] "6 60 4.59999999999999"
# [1] "7 70 4.84"
# [1] "8 80 5.07999999999998"
# [1] "9 90 5.29999999999998"
# [1] "10 100 5.52000000000001"
# [1] "11 110 5.72"
# [1] "12 120 5.97"
# [1] "13 130 6.27000000000001"


# [1] "1 30 4.51999999999998"
# [1] "2 35 4.41999999999999"
# [1] "3 40 4.45000000000002"
# [1] "4 45 4.44999999999999"
# [1] "5 50 4.49000000000001"
# [1] "6 55 4.59"



# With 4D data

# [1] "1 1 5.66"
# [1] "2 5 1.13"
# [1] "3 10 0.670000000000016"
# [1] "4 20 0.460000000000008"
# [1] "5 30 0.430000000000007"
# [1] "6 40 0.689999999999998"
# [1] "7 50 0.490000000000009"
# [1] "8 60 0.450000000000017"
# [1] "9 70 0.560000000000002"
# [1] "10 80 0.52000000000001"
# [1] "11 90 0.5"
# [1] "12 100 0.580000000000013"
# [1] "13 110 0.550000000000011"
# [1] "14 120 0.579999999999984"
# [1] "15 130 0.719999999999999"
# [1] "16 300 1.08000000000001"
# [1] "17 1000 3.19999999999999"





# Check of split_speed, 50x speedup here
gp <- GauPro(x,y, verbose=2);gp$theta
# Make sure predictions are the same
summary(gp$pred(XX, se.fit=T, split_speed=F) - gp$pred(XX, se.fit=T, split_speed=T))
# 30x speedup on full XX
system.time(gp$pred(XX, se.fit=T, split_speed=T))
system.time(gp$pred(XX, se.fit=T, split_speed=F))

# Slower on 160 by a little
microbenchmark::microbenchmark(gp$pred(XX[1:161,], se.fit=T, split_speed=F),gp$pred(XX[1:161,], se.fit=T, split_speed=T))
microbenchmark::microbenchmark(gp$pred(XX[1:160,], se.fit=T, split_speed=F),gp$pred(XX[1:160,], se.fit=T, split_speed=T))
# 5x faster on 1610
microbenchmark::microbenchmark(gp$pred(XX[1:1610,], se.fit=T, split_speed=F),gp$pred(XX[1:1610,], se.fit=T, split_speed=T))
# About the same at 200
microbenchmark::microbenchmark(gp$pred(XX[1:200,], se.fit=T, split_speed=F),gp$pred(XX[1:200,], se.fit=T, split_speed=T))

