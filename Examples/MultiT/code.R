library(LaplacesDemon)
library(pracma)
library(SQUAREM)
library(quasiNewtonMM)
source("qnamm.r")

VecToMat <- function(vec, dim){
  a <- matrix(0, dim, dim)
  a[upper.tri(a, diag = TRUE)] <- vec
  a = a + t(a)
  diag(a) <- diag(a)/2
  return(a)
}

likelihood <- function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  like <- (n*log(det(sigma)))/2
  sig.inv <- solve(sigma)
  for (i in 1:n){
    like = like - ((1+dim)/2) * log(1 + t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
  }
  return(like)
}

update <-  function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  sig.inv <- solve(sigma)
  weights <- as.matrix(rep(0, n))

  new.mu <- rep(0,dim)
  new.sigma <- matrix(0, dim, dim)

  for (i in 1:n){
    weights[i] <- (1+dim)/(1 + (t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu)))
    new.mu <- new.mu + weights[i]*as.matrix(data[i,])
    new.sigma <- new.sigma + weights[i]*((data[i,] - mu) %*% t(data[i,] - mu))
  }
  new.mu <- new.mu/sum(weights)
  new.sigma <- new.sigma/n

  return (c(new.mu, upper.triangle(new.sigma, diag = TRUE)))
}

update_pxem <-  function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  sig.inv <- solve(sigma)
  weights <- as.matrix(rep(0, n))

  new.mu <- rep(0,dim)
  new.sigma <- matrix(0, dim, dim)

  for (i in 1:n){
    weights[i] <- (1+dim)/(1 + (t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu)))
    new.mu <- new.mu + weights[i]*as.matrix(data[i,])
    new.sigma <- new.sigma + weights[i]*((data[i,] - mu) %*% t(data[i,] - mu))
  }
  new.mu <- new.mu/sum(weights)
  new.sigma <- new.sigma/sum(weights)

  return (c(new.mu, upper.triangle(new.sigma, diag = TRUE)))
}


##################################################

set.seed(10)
dim <- 10
tol <- 1e-7
P <- (dim/2)*(dim+3)
n <- 100
mu <- rep(0, dim)
u <- matrix(rnorm(dim*dim), dim, dim)
sigma <- t(u) %*% u
N <- 10
start_rep <- matrix(0, nrow = N, ncol = P)
for (i in 1:N)
{
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start_rep[i,] <- c(mu0, upper.triangle(sigma0, diag=  TRUE))
}


###########################################
## EM Algorithm
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
    if(iter %% 1000 == 0) print(iter)
    new <- update(now, n=n, dim=dim, data=data)
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- likelihood(new, n=n, dim=dim, data=data)
  eval_mm[i] <- iter

}

print(quantile(time_mm, probs = c(0, .5, 1)))
print(quantile(eval_mm, probs = c(0, .5, 1)))
print(quantile(obj_mm, probs = c(0, .5, 1)))

###########################################
## PX-EM
###########################################

time_pxem <- rep(0, N)
obj_pxem <- rep(0, N)
eval_pxem <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]

  start.time <- Sys.time()
  fp <- fpiter(par = start, n=n, dim=dim, data=data, fixptfn = update_pxem,
               objfn = likelihood, control = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()

  if(fp$convergence){
    time_pxem[i] <- end.time - start.time
    obj_pxem[i] <- fp$value.objfn
    eval_pxem[i] <- fp$fpevals
  } else{
    time_pxem[i] <- NA
    obj_pxem[i] <- NA
    eval_pxem[i] <- NA
  }
}

print(quantile(time_pxem, probs = c(0, .5, 1)))
print(quantile(eval_pxem, probs = c(0, .5, 1)))
print(quantile(obj_pxem, probs = c(0, .5, 1)))

##########################################
## ZAL, q=1
##########################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = likelihood, qn=2, n=n, dim=dim, data=data, max_iter = 1e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
}

print(quantile(time_zal, probs = c(0, .5, 1)))
print(quantile(eval_zal, probs = c(0, .5, 1)))
print(quantile(obj_zal, probs = c(0, .5, 1)))

########################################
## BQN1
########################################

time_bqn1 <- rep(0, N)
obj_bqn1 <- rep(0, N)
eval_bqn1 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
             control = list(qn=1, tol = tol, objfn.inc = 1, step.max = 1e7, maxiter = 5e4))
  end.time <- Sys.time()

  time_bqn1[i] <- end.time - start.time
  obj_bqn1[i] <- fp$value.objfn
  eval_bqn1[i] <- fp$fpevals
}

print(quantile(time_bqn1, probs = c(0, .5, 1)))
print(quantile(eval_bqn1, probs = c(0, .5, 1)))
print(quantile(obj_bqn1, probs = c(0, .5, 1)))

########################################
## BQN2
########################################

time_bqn2 <- rep(0, N)
obj_bqn2 <- rep(0, N)
eval_bqn2 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
             control = list(qn=2, tol = tol, objfn.inc = 1, step.max = 1e6, maxiter = 5e4))
  end.time <- Sys.time()
  
  time_bqn2[i] <- end.time - start.time
  obj_bqn2[i] <- fp$value.objfn
  eval_bqn2[i] <- fp$fpevals
}

print(quantile(time_bqn2, probs = c(0, .5, 1)))
print(quantile(eval_bqn2, probs = c(0, .5, 1)))
print(quantile(obj_bqn2, probs = c(0, .5, 1)))

########################################
## L-BQN
########################################

time_lbqn <- rep(0, N)
obj_lbqn <- rep(0, N)
eval_lbqn <- rep(0, N)


for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
              control = list(m=10, tol = tol, objfn.inc = 1, maxiter = 5e4))
  end.time <- Sys.time()
  time_lbqn[i] <- end.time - start.time
  obj_lbqn[i] <- fp$value.objfn
  eval_lbqn[i] <- fp$fpevals
  }
print(quantile(time_lbqn, probs = c(0, .5, 1)))
print(quantile(eval_lbqn, probs = c(0, .5, 1)))
print(quantile(obj_lbqn, probs = c(0, .5, 1)))


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
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
}

print(quantile(time_sq1, probs = c(0, .5, 1)))
print(quantile(eval_sq1, probs = c(0, .5, 1)))
print(quantile(obj_sq1, probs = c(0, .5, 1)))

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
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
}
print(quantile(time_sq2, probs = c(0, .5, 1)))
print(quantile(eval_sq2, probs = c(0, .5, 1)))
print(quantile(obj_sq2, probs = c(0, .5, 1)))

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
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood,n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
}

print(quantile(time_sq3, probs = c(0, .5, 1)))
print(quantile(eval_sq3, probs = c(0, .5, 1)))
print(quantile(obj_sq3, probs = c(0, .5, 1)))



save(time_mm, time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bqn1, eval_bqn2, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bqn1, obj_bqn2, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = "Out/multiT-objects.Rdata")


