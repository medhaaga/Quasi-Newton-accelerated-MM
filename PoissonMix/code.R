set.seed(100)
library(pracma)
library(BfgsQN)
library(SQUAREM)
source("qnamm.r")

likelihood <- function(par, data)
{
  p <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]
  like <- 0
  for (i in 0:9){
    like <- like + data[i+1,2] * log(((p * exp(-mu1) * (mu1 ^ i))/factorial(i)) + (((1-p) * exp(-mu2) * (mu2 ^ i))/factorial(i)))
  }
  return(-like)
}


update <- function(par, data){

  p <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]

  pi <- matrix(0, nrow = 10, ncol = 2)
  denom <- rep(0, 10)
  for (i in 0:9){
    denom[i+1] <- (p * (mu1 ^ i) * exp(-mu1)) + ((1-p) * (mu2 ^ i) * exp(-mu2))
    pi[i+1,1] <- (p * (mu1 ^ i) * exp(-mu1))/denom[i+1]
    pi[i+1,2] <- 1 - pi[i+1,1]
  }

  new_p <- dot(data[,2], pi[,1])/sum(data[,2])
  new_mu1 <- dot((data[,1]*data[,2]), pi[,1])/dot(data[,2], pi[,1])
  new_mu2 <- dot((data[,1]*data[,2]), pi[,2])/dot(data[,2], pi[,2])
  return(c(new_p, new_mu1, new_mu2))
}

############################################################3

set.seed(10)
N <- 1
deaths <- (0:9)
freq <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
data <- as.matrix(cbind(deaths, freq))
truth <- c(.3599, 1.256, 2.663)
#start_rep <- matrix(c(round(runif(10, .1, .9), 2), round(runif(10, 1, 50), 2), round(runif(10, 1, 50), 2)), nrow = N, ncol = 3)
start_rep <- matrix(c(.2870, 1.101, 2.582), nrow = N)
dim <- 3
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
  chain <- matrix(0, nrow = 1e5, ncol = dim)

  while((diff > tol))
  {
    iter <- iter + 1
    if(iter %% 2 == 0) print(iter)
    new <- update(now, data)
    chain[iter,] <- new
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- likelihood(new, data)
  eval_mm[i] <- iter
  chain <- chain[1:iter,]
  print(chain[iter,])
}

print(quantile(time_mm, c(0.25, .5, 0.75)))
print(quantile(obj_mm, c(0.25, .5, 0.75)))
print(quantile(eval_mm, c(0.25, .5, 0.75)))

########################################
## Classical BFGS, q=1
########################################

time_bfgs1 <- rep(0, N)
obj_bfgs1 <- rep(0, N)
eval_bfgs1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BFGS(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(qn=1, tol = tol, objfn.inc = 1, step.max = 10, maxiter = 5e3, intermed = TRUE))
  end.time <- Sys.time()

  time_bfgs1[i] <- end.time - start.time
  obj_bfgs1[i] <- fp$value.objfn
  eval_bfgs1[i] <- fp$fpevals
  chain <- (fp$p.inter)
  print(chain[fp$iter,])
}

print(quantile(time_bfgs1, c(0.25, .5, 0.75)))
print(quantile(obj_bfgs1, c(0.25, .5, 0.75)))
print(quantile(eval_bfgs1, c(0.25, .5, 0.75)))



########################################
## Classical BFGS, q=2
########################################

time_bfgs2 <- rep(0, N)
obj_bfgs2 <- rep(0, N)
eval_bfgs2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BFGS(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(qn=2, tol = tol, objfn.inc = 1, step.max = 10, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_bfgs2[i] <- end.time - start.time
  obj_bfgs2[i] <- fp$value.objfn
  eval_bfgs2[i] <- fp$fpevals
  print(chain[fp$iter,])
}

print(quantile(time_bfgs2, c(0.25, .5, 0.75)))
print(quantile(obj_bfgs2, c(0.25, .5, 0.75)))
print(quantile(eval_bfgs2, c(0.25, .5, 0.75)))

########################################
## L-BFGS
########################################

time_lbfgs <- rep(0, N)
obj_lbfgs <- rep(0, N)
eval_lbfgs <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBFGS(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(m=10, tol = tol, objfn.inc = 1, maxiter = 1e4, intermed = TRUE))
  end.time <- Sys.time()
  time_lbfgs[i] <- end.time - start.time
  obj_lbfgs[i] <- fp$value.objfn
  eval_lbfgs[i] <- fp$fpevals
}

print(quantile(time_lbfgs, c(0, .5, 1)))
print(quantile(obj_lbfgs, c(0, .5, 1)))
print(quantile(eval_lbfgs, c(0, .5, 1)))


##########################################
### SqS1
#############################################

time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)
eval_sq1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
  print(fp$convergence)
}

print(quantile(time_sq1, c(0, .5, 1)))
print(quantile(obj_sq1, c(0, .5, 1)))
print(quantile(eval_sq1, c(0, .5, 1)))

##########################################
### SqS2
#############################################

time_sq2 <- rep(0, N)
obj_sq2 <- rep(0, N)
eval_sq2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
}

print(quantile(time_sq2, c(0, .5, 1)))
print(quantile(obj_sq2, c(0, .5, 1)))
print(quantile(eval_sq2, c(0, .5, 1)))


##########################################
### SqS3
#############################################

time_sq3 <- rep(0, N)
obj_sq3 <- rep(0, N)
eval_sq3 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, data=data, control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
}

print(quantile(time_sq3, c(0, .5, 1)))
print(quantile(obj_sq3, c(0, .5, 1)))
print(quantile(eval_sq3, c(0, .5, 1)))

##########################################
## Zhou's quasi-Newton for q=1
##########################################

time_zal <- rep(0, N)
obj_zal <- rep(0, N)
eval_zal <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = likelihood, qn=3, data=data, max_iter = 1e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
}

print(quantile(time_zal, c(0, .5, 1)))
print(quantile(obj_zal, c(0, .5, 1)))
print(quantile(eval_zal, c(0, .5, 1)))

