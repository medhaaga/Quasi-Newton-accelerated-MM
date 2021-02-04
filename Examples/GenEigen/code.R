library(pracma)
library(SQUAREM)
library(quasiNewtonMM)
source("qnamm.r")

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

N <- 1
start_rep <- matrix(rnorm(N*dim, mean = 1, sd = 10), nrow = N, ncol = dim)
dir <- "descent"
tol <- 1e-7


###########################################
## Unaccelerated MM Algorithm
###########################################

time_mm <- rep(0, N)
obj_mm <- rep(0, N)
eval_mm <- rep(0, N)

for (i in 1:N){
  print(i)
  now <- start_rep[i,]
  new <- start_rep[i,]
  diff <- 100

  iter <- 0
  start.time <- Sys.time()
  while((diff > tol))
  {
    iter <- iter + 1
    new <- update(now, A, B, dir = dir)
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- rayleigh(new, A, B, dir)
  eval_mm[i] <- iter
}

print(quantile(time_mm, probs = c(.1, .5, .9)))
print(quantile(eval_mm, probs = c(.1, .5, .9)))
print(quantile(obj_mm, probs = c(.1, .5, .9)))


########################################
## BQN, q=1
########################################

time_bqn1 <- rep(0, N)
obj_bqn1 <- rep(0, N)
eval_bqn1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(qn=1, tol = tol, objfn.inc = 1, step.max = 1000, step.min=.1, maxiter = 5e4))
  end.time <- Sys.time()

  time_bqn1[i] <- end.time - start.time
  obj_bqn1[i] <- fp$value.objfn
  eval_bqn1[i] <- fp$fpevals
}

print(quantile(time_bqn1, probs = c(.1, .5, .9)))
print(quantile(eval_bqn1, probs = c(.1, .5, .9)))
print(quantile(obj_bqn1, probs = c(.1, .5, .9)))

########################################
## BQN, q=2
########################################


time_bqn2 <- rep(0, N)
obj_bqn2 <- rep(0, N)
eval_bqn2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(qn=2, tol = tol, objfn.inc = 1, step.max = 1000, step.min=1, maxiter = 5e4))
  end.time <- Sys.time()
  
  time_bqn2[i] <- end.time - start.time
  obj_bqn2[i] <- fp$value.objfn
  eval_bqn2[i] <- fp$fpevals
}


print(quantile(time_bqn2, probs = c(.1, .5, .9)))
print(quantile(eval_bqn2, probs = c(.1, .5, .9)))
print(quantile(obj_bqn2, probs = c(.1, .5, .9)))

########################################
## L-BFGS
########################################

time_lbqn <- rep(0, N)
obj_lbqn <- rep(0, N)
eval_lbqn <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(m=10, tol = tol, objfn.inc = 1, maxiter = 5e4))
  end.time <- Sys.time()
  time_lbqn[i] <- end.time - start.time
  obj_lbqn[i] <- fp$value.objfn
  eval_lbqn[i] <- fp$fpevals

}
print(quantile(time_lbqn, probs = c(.1, .5, .9)))
print(quantile(eval_lbqn, probs = c(.1, .5, .9)))
print(quantile(obj_lbqn, probs = c(.1, .5, .9)))


##########################################
### SqS1
#############################################

time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)
eval_sq1 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 1, maxiter = 5e4))
  end.time <- Sys.time()

  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
  }

print(quantile(time_sq1, probs = c(.1, .5, .9)))
print(quantile(eval_sq1, probs = c(.1, .5, .9)))
print(quantile(obj_sq1, probs = c(.1, .5, .9)))

##########################################
### SqS2
#############################################

time_sq2 <- rep(0, N)
obj_sq2 <- rep(0, N)
eval_sq2 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 2, maxiter = 5e4))
  end.time <- Sys.time()

  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
  }

print(quantile(time_sq2, probs = c(.1, .5, .9)))
print(quantile(eval_sq2, probs = c(.1, .5, .9)))
print(quantile(obj_sq2, probs = c(.1, .5, .9)))

##########################################
### SqS3
#############################################

time_sq3 <- rep(0, N)
obj_sq3 <- rep(0, N)
eval_sq3 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K=2, tol = tol, method = 3, maxiter = 5e4))
  end.time <- Sys.time()

  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
  }

print(quantile(time_sq3, probs = c(.1, .5, .9)))
print(quantile(eval_sq3, probs = c(.1, .5, .9)))
print(quantile(obj_sq3, probs = c(.1, .5, .9)))

##########################################
## Zhou's quasi-Newton for q=2
##########################################

time_zal <- rep(0, N)
obj_zal <- rep(0, N)
eval_zal <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = rayleigh, qn=2, A=A, B=B, dir=dir, max_iter = 5e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
  
}

print(quantile(time_zal, probs = c(.1, .5, .9)))
print(quantile(eval_zal, probs = c(.1, .5, .9)))
print(quantile(obj_zal, probs = c(.1, .5, .9)))

##############################################

save(time_mm, time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bqn1, eval_bqn2, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bqn1, obj_bqn2, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = "Out/eigen-objects_sd10.Rdata")

load("Out/eigen-objects_sd10.Rdata")
time_range <- range(time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal)
eval_range <- range(eval_bqn1, eval_bqn2, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal)
obj_range <- range(obj_bqn1, obj_bqn2, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal)

plot(time_bqn1, eval_bqn1, pch=19, col  ="red", xlim = time_range, ylim = eval_range)
points(time_bqn2, eval_bqn2, pch=19, col = "green3")
points(time_lbqn, eval_lbqn, pch=19, col = "blue")
points(time_sq1, eval_sq1, pch=19, col = "pink")
points(time_sq2, eval_sq2, pch=19)
points(time_sq3, eval_sq3, pch=19)
points(time_zal, eval_zal, pch=19)
