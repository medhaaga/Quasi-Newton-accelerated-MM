##################################################
############### Probit Regression ################
##################################################

set.seed(1)
rm(list = ls())
library(pracma)
library(quasiNewtonMM)
library(turboEM)
library(daarem)
library(RColorBrewer)
library(SQUAREM)

update <- function(b, X, Y, foo)
{
  X_beta <- X %*% b
  U <- matrix(0, n, 1)
  U[Y==1] <-   X_beta[Y==1] + dnorm(-X_beta[Y==1])/(1 - pnorm(-X_beta[Y==1]))
  U[Y==0] <-   X_beta[Y==0] - dnorm(X_beta[Y==0])/(1 - pnorm(X_beta[Y==0]))
  return(foo %*% t(X) %*% U)
}

neg_loglikelihood <- function(b, X, Y, foo)
{
  X_beta <- X %*% b
  phi_X_beta <- pnorm(X_beta)
  like <- prod((phi_X_beta^Y)*((1-phi_X_beta)^(1-Y)))
  return(-log(like))
}

loglikelihood <- function(b, X, Y, foo)
{
  X_beta <- X %*% b
  phi_X_beta <- pnorm(X_beta)
  like <- prod((phi_X_beta^Y)*((1-phi_X_beta)^(1-Y)))
  return(log(like))
}

n <- 2000
p <- 1000
X <- matrix(rnorm(n*p), n, p)
true_beta <- matrix(rt(p, df=2)/2 + 2, p, 1)
Z <- X %*% true_beta + matrix(rnorm(n), n, 1)
Y <- matrix(0, n, 1)
Y[Z>0] <- 1
foo <- solve(t(X) %*% X)
true_neg_loglikelihood <- neg_loglikelihood(true_beta, X, Y)


N <- 1
start.all <- 0.1*matrix(runif(N*p), nrow = N, ncol = p)
tol = 1e-7

time_mm <- rep(0,N)
eval_mm <- rep(0,N)
obj_mm <- rep(0,N)

for (j in 1:N)
{
  print(j)

  start <- as.matrix(start.all[j,])
  now <- start
  new <- start
  iter <- 1
  diff <- 100


  start.time <- Sys.time()
  while(diff > tol)
  {
    new <- update(now, X, Y, foo)
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
    if (iter %% 100 == 0){
      print(neg_loglikelihood(new, X, Y, foo))
    }
  }
  end.time <- Sys.time()
  obj_mm[j] <- neg_loglikelihood(new, X, Y, foo)
  time_mm [j] <- end.time - start.time
  eval_mm[j] <- iter
}

print(round(quantile(time_mm, c(.5, .25, .75)), 3))
print(quantile(eval_mm, c(.5, .25, .75)))
print(quantile(obj_mm, c(.5, .25, .75)))

#####################################################
#### L-BQN
#####################################################


time_lbqn <- rep(0, N)
obj_lbqn <- rep(0, N)
eval_lbqn <- rep(0, N)


for (i in 1:N)
{
  print(i)
  start <-as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- LBQN(par = start, X=X, Y=Y, foo=foo, fixptfn = update, objfn = neg_loglikelihood,
             control = list(tol = 1e-7, maxiter = 1e5, m = min(p,10), verbose=TRUE))
  end.time <- Sys.time()
  
  
  time_lbqn[i] <- end.time - start.time
  obj_lbqn[i] <- fp$value.objfn
  eval_lbqn[i] <- fp$fpevals
}

print(round(quantile(time_lbqn, c(.5, .25, .75)), 3))
print(quantile(eval_lbqn, c(.5, .25, .75)))
print(round(quantile(obj_lbqn, c(.5, .25, .75)), 6))

##########################################
#### SQUAREM - 1
##########################################

eval_sq1 <- rep(0, N)
time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = neg_loglikelihood, X=X, Y=Y, foo=foo, method = "squarem", control.method = list(K=1, version=1), control.run = list(tol=tol, maxiter=1e5))
  end.time <- Sys.time()
  
  
  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpeval
}

print(round(quantile(time_sq1, c(.5, .25, .75)), 3))
print(quantile(eval_sq1, c(.5, .25, .75)))
print(round(quantile(obj_sq1, c(.5, .25, .75)), 6))

#################################################
##### ZAL; q=2
#################################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = neg_loglikelihood, X=X, Y=Y, foo=foo, method = "qn", control.method = list(qn=2), control.run = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  print(fp$convergence)
  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$value.objfn
  eval_zal[i] <- fp$fpeval
}

print(round(quantile(time_zal, c(.5, .25, .75)), 3))
print(quantile(eval_zal, c(.5, .25, .75)))
print(round(quantile(obj_zal, c(.5, .25, .75)), 6))


#################################################
##### ZAL; q=min(p,10)
#################################################

time_zal2 <- rep(NA, N)
obj_zal2 <- rep(NA, N)
eval_zal2 <- rep(NA, N)

for (i in 1:N){
  start <- as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = neg_loglikelihood, X=X, Y=Y, foo=foo, method = "qn", control.method = list(qn=min(p,10)), control.run = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  print(fp$convergence)
  time_zal2[i] <- end.time - start.time
  obj_zal2[i] <- fp$value.objfn
  eval_zal2[i] <- fp$fpeval
}

print(round(quantile(time_zal2, c(.5, .25, .75)), 3))
print(quantile(eval_zal2, c(.5, .25, .75)))
print(round(quantile(obj_zal2, c(.5, .25, .75)), 6))


#################################################
##### DAAREM
#################################################

time_dar <- rep(NA, N)
obj_dar <- rep(NA, N)
eval_dar <- rep(NA, N)

for (i in 1:N){
  print(paste("Iter: ", i))
  start <- as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- daarem(par = start, fixptfn = update, objfn = loglikelihood, X=X, Y=Y, foo=foo, control = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  time_dar[i] <- end.time - start.time
  obj_dar[i] <- neg_loglikelihood(fp$par, X, Y, foo)
  eval_dar[i] <- fp$fpeval
}

print(round(quantile(time_dar, c(.5, .25, .75)), 3))
print(quantile(eval_dar, c(.5, .25, .75)))
print(round(quantile(obj_dar, c(.5, .25, .75)), 7))

