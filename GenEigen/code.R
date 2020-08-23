library(pracma)
library(SQUAREM)
library(BfgsQN)
source("qnamm.r")
library(ggfortify)

rayleigh <- function(x, A, B, dir){
  x <- as.matrix(x)
  num <- t(x) %*% A %*% x
  denom <- t(x) %*% B %*% x
  if (dir == "descent")
    return(num/denom) else
      return(-num/denom)
}

update <- function(x, A, B, dir = c("ascent", "descent")){
  x <- as.matrix(x)
  u <- as.matrix(x)
  v <- (A - as.numeric((t(x) %*% A %*% x)/(t(x) %*% B %*% x))*B) %*% x

  uAu <- t(u) %*% A %*% u
  vAv <- t(v) %*% A %*% v
  uAv <- t(u) %*% A %*% v
  vAu <- t(v) %*% A %*% u
  uBu <- t(u) %*% B %*% u
  vBv <- t(v) %*% B %*% v
  uBv <- t(u) %*% B %*% v
  vBu <- t(v) %*% B %*% u

  a <- (vAv*vBu + vAu*vBv + vAv*uBv) - (vAv*vBu + uAv*vBv + vAu*vBv)
  b <- (vAv*uBu + vAu*vBu + vAu*uBv) - (vAu*vBu + uAv*vBu + uAu*vBv)
  c <- (vAu*uBu) - (uAu*vBu)

  delta <- (b^2 - 4*a*c)
  if(delta > 0){ # first case D>0
    x_1 = (-b+sqrt(delta))/(2*a)
    x_2 = (-b-sqrt(delta))/(2*a)
    C <- c(x_1,x_2)
  }
  else if(delta == 0){ # second case D=0
    x = -b/(2*a)
    C <- (c(x,x))
  }
  else {stop("No roots")} # third case D<0

  x1 <- u + C[1]*v
  x2 <- u + C[2]*v
  num <- t(x1) %*% A %*% x1
  denom <- t(x1) %*% B %*% x1
  R1 <- num/denom
  num <- t(x2) %*% A %*% x2
  denom <- t(x2) %*% B %*% x2
  R2 <- num/denom

  if(R1 > R2) {
    ascent <- x1
    descent <- x2
  }
  else {
    ascent <- x2
    descent <- x1
  }

  if (dir == "ascent") {return(ascent)}
  else {return(descent)}

}

set.seed(1)
dim <- 100
C <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
D <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
A <- C + t(C)
B <- D %*% t(D)

N <- 10
start_rep <- matrix(rnorm(N*dim, mean = 1, sd = 10), nrow = N, ncol = dim)
dir <- "descent"
tol <- 1e-7



###########################################
## Naive Algorithm
###########################################

time_mm <- rep(0, N)
obj_mm <- rep(0, N)
eval_mm <- rep(0, N)

for (i in 1:N){
  now <- start_rep[i,]
  new <- start_rep[i,]
  diff <- 100

  iter <- 0
  start.time <- Sys.time()
  chain <- matrix(0, nrow = 1e6, ncol = dim)
  while((diff > tol))
  {
    iter <- iter + 1
    if(iter %% 1000 == 0) print(iter)
    new <- update(now, A, B, dir = dir)
    chain[iter,] <- new
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- rayleigh(new, A, B, dir)
  eval_mm[i] <- iter
  chain <- chain[1:iter,]

  plot(chain[,27], chain[,93], ylim=range(chain[,93]), xlim = range(chain[,27]),
       xlab="PC 1", ylab="PC 2")
  points(chain[iter,27], chain[iter, 93], col = "red", pch = 19, cex=2)
  points(chain[iter,27], chain[1, 93], col = "green", pch = 19, cex=2)
}

########################################
## Classical BFGS
########################################


time_bfgs <- rep(0, N)
obj_bfgs <- rep(0, N)
eval_bfgs <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BFGS(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(tol = tol, objfn.inc = 1e-2, obj.tol = 1e-9, step.max = 10, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_bfgs[i] <- end.time - start.time
  obj_bfgs[i] <- fp$value.objfn
  eval_bfgs[i] <- fp$fpevals

  chain <- (fp$p.inter)
  plot(chain[,27], chain[,93], ylim=range(chain[,93]), xlim = range(chain[,27]),
       xlab="PC 1", ylab="PC 2")
  points(chain[dim(fp$p.inter)[1],27], chain[dim(fp$p.inter)[1], 93], col = "red", pch = 19, cex=2)
  points(chain[1,27], chain[1, 93], col = "green", pch = 19, cex=2)
}

########################################
## L-BFGS
########################################

time_lbfgs <- rep(0, N)
obj_lbfgs <- rep(0, N)
eval_lbfgs <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBFGS(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(m=10, tol = tol, objfn.inc = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  time_lbfgs[i] <- end.time - start.time
  obj_lbfgs[i] <- fp$value.objfn
  eval_lbfgs[i] <- fp$fpevals

  chain <- (fp$p.inter)
  plot(chain[,27], chain[,93], ylim=range(chain[,93]), xlim = range(chain[,27]),
       xlab="PC 1", ylab="PC 2")
  points(chain[dim(fp$p.inter)[1],27], chain[dim(fp$p.inter)[1], 93], col = "red", pch = 19, cex=2)
  points(chain[1,27], chain[1, 93], col = "green", pch = 19, cex=2)
}


##########################################
### SqS1
#############################################

time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)
eval_sq1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
  }

##########################################
### SqS2
#############################################

time_sq2 <- rep(0, N)
obj_sq2 <- rep(0, N)
eval_sq2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
  }

##########################################
### SqS3
#############################################

time_sq3 <- rep(0, N)
obj_sq3 <- rep(0, N)
eval_sq3 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
  }

##########################################
## Zhou's quasi-Newton for q=1
##########################################

time_zal <- rep(0, N)
obj_zal <- rep(0, N)
eval_zal <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = rayleigh, qn=2, A=A, B=B, dir=dir, max_iter = 1e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals

  chain <- t(fp$Xhist)
  plot(chain[,27], chain[,93], ylim=range(chain[,93]), xlim = range(chain[,27]),
       xlab="PC 1", ylab="PC 2")
  points(chain[dim(chain)[1],27], chain[dim(chain)[1], 93], col = "red", pch = 19, cex=2)
  points(chain[1,27], chain[1, 93], col = "green", pch = 19, cex=2)
}

##############################################

save(time_mm, time_bfgs, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bfgs, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bfgs, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = "Out/objects_sd10.Rdata")

time_range <- range(time_mm, time_bfgs, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal)
eval_range <- range(eval_mm, eval_bfgs, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal)
obj_range <- range(obj_mm, obj_bfgs, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal)

pdf(file = "Out/eigen-objVSeval_sd10.pdf")
plot(eval_mm, obj_mm, col = "red", xlim=eval_range, ylim = obj_range, pch=19, cex=2, xlab="Evaluations", ylab="Objective")
points(eval_bfgs, obj_bfgs, col="purple", pch=19, cex=2)
points(eval_lbfgs, obj_lbfgs, col="pink", pch=19, cex=2)
points(eval_sq1, obj_sq1, col="lightblue", pch=19, cex=2)
#points(eval_sq2, obj_sq2, col=6, pch=19)
#points(eval_sq3, obj_sq3, col=7, pch=19)
points(eval_zal, obj_zal, col="steelblue1", pch=19, cex=2)
legend("topright", legend = c("MM", "BFGS", "L-BFGS", "SQUAREM", "ZAL"),
       col =c("red", "purple", "pink", "lightblue", "steelblue1"), pch=19, cex=2)
dev.off()

pdf(file = "Out/eigen-objVStime_sd10.pdf")
plot(time_mm, obj_mm, col = "red", xlim=time_range, ylim = obj_range, pch=19, cex=2, xlab="Time", ylab="Objective")
points(time_bfgs, obj_bfgs, col="purple", pch=19, cex=2)
points(time_lbfgs, obj_lbfgs, col="pink", pch=19, cex=2)
points(time_sq1, obj_sq1, col="lightblue", pch=19, cex=2)
#points(eval_sq2, obj_sq2, col=6, pch=19)
#points(eval_sq3, obj_sq3, col=7, pch=19)
points(time_zal, obj_zal, col="steelblue1", pch=19, cex=2)
legend("topright", legend = c("MM", "BFGS", "L-BFGS", "SQUAREM", "ZAL"),
       col =c("red", "purple", "pink", "lightblue", "steelblue1"), pch=19, cex=2)
dev.off()

