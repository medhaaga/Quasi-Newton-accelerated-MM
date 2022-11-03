
#############################################
######## Multivariate t-distribution ########
#############################################

rm(list = ls())
library(LaplacesDemon)
library(pracma)
library(turboEM)
library(quasiNewtonMM)
source("functions.R")
##################################################

set.seed(1)
dim <- 25
tol <- 1e-7
P <- (dim/2)*(dim+3)
n <- 1000
mu <- rep(0, dim)
u <- matrix(rnorm(dim*dim, sd = 1), dim, dim)
sigma <- t(u) %*% u
N <- 10
start_rep <- matrix(0, nrow = N, ncol = P)
data <- array(0, dim = c(N, n, dim))

for (i in 1:N)
{
  data[i,,] <- (rmvc(n=n, mu = mu, S = sigma))
  mu0 <- colMeans(data[i,,])
  sigma0 <- cov(data[i,,])
  start_rep[i,] <- c(mu0, upper.triangle(sigma0, diag=  TRUE))
}


###########################################
## Unaccelerated EM Algorithm
###########################################

time_mm <- rep(NA, N)
obj_mm <- rep(NA, N)
eval_mm <- rep(NA, N)

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
    if(iter %% 1000 == 0) print(diff)
    new <- update(now, n=n, dim=dim, data=data[i,,])
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- likelihood(new, n=n, dim=dim, data=data[i,,])
  eval_mm[i] <- iter

}

print(paste("Number of failures:", sum(is.na(time_mm))))
print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

###########################################
## PX-EM
###########################################

time_pxem <- rep(NA, N)
obj_pxem <- rep(NA, N)
eval_pxem <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]

  start.time <- Sys.time()
  fp <- fpiter(par = start, n=n, dim=dim, data=data[i,,], fixptfn = update_pxem,
               objfn = likelihood, control = list(tol = tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_pxem[i] <- end.time - start.time
    obj_pxem[i] <- fp$value.objfn
    eval_pxem[i] <- fp$fpevals
  }
}

print(paste("Number of failures:", sum(is.na(time_pxem))))
print(round(quantile(time_pxem[!is.na(time_pxem)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_pxem[!is.na(eval_pxem)], probs = c(.5, .25, .75)))
print(round(quantile(obj_pxem[!is.na(obj_pxem)], probs = c(.5, .25, .75)), 4))


########################################
## BQN, q=1
########################################

time_bqn1 <- rep(NA, N)
obj_bqn1 <- rep(NA, N)
eval_bqn1 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], control = list(qn=1, tol = tol, maxiter = 5e4, objfn.inc = .1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_bqn1[i] <- end.time - start.time
    obj_bqn1[i] <- fp$value.objfn
    eval_bqn1[i] <- fp$fpevals
  }
}

print(paste("Number of failures:", sum(is.na(time_bqn1))))
print(round(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)), 4))

########################################
## BQN, q=2
########################################

time_bqn2 <- rep(NA, N)
obj_bqn2 <- rep(NA, N)
eval_bqn2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,],
             control = list(qn=2, tol = tol, maxiter = 5e3, objfn.inc = .1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_bqn2[i] <- end.time - start.time
    obj_bqn2[i] <- fp$value.objfn
    eval_bqn2[i] <- fp$fpevals
  }
}

print(paste("Number of failures:", sum(is.na(time_bqn2))))
print(round(quantile(time_bqn2[!is.na(time_bqn2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn2[!is.na(eval_bqn2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn2[!is.na(obj_bqn2)], probs = c(.5, .25, .75)), 4))

########################################
## BQN, q=min(p,10)
########################################

time_bqn3 <- rep(NA, N)
obj_bqn3 <- rep(NA, N)
eval_bqn3 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- try(BQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,],
            control = list(qn=min(P,3), tol = tol, maxiter = 5e3, objfn.inc = .1)))
  end.time <- Sys.time()
  if(inherits(fp, "try-error")){
  print("Outside param space")
  }
  else if(fp$convergence){
    time_bqn3[i] <- end.time - start.time
    obj_bqn3[i] <- fp$value.objfn
    eval_bqn3[i] <- fp$fpevals
  }
}

print(paste("Number of failures:", sum(is.na(time_bqn3))))
print(round(quantile(time_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn3[!is.na(eval_bqn3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn3[!is.na(obj_bqn3)], probs = c(.5, .25, .75)), 4))

########################################
## L-BQN
########################################

time_lbqn <- rep(NA, N)
obj_lbqn <- rep(NA, N)
eval_lbqn <- rep(NA, N)


for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBQN(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,],
              control = list(m=min(P, 10), tol = tol, maxiter = 1e3))
  end.time <- Sys.time()
  if(fp$convergence){
    time_lbqn[i] <- end.time - start.time
    obj_lbqn[i] <- fp$value.objfn
    eval_lbqn[i] <- fp$fpevals
  }
}

print(paste("Number of failures:", sum(is.na(time_lbqn))))
print(round(quantile(time_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_lbqn[!is.na(eval_lbqn)], probs = c(.5, .25, .75)))
print(round(quantile(obj_lbqn[!is.na(obj_lbqn)], probs = c(.5, .25, .75)), 4))


##########################################
### SqS1
#############################################

time_sq1 <- rep(NA, N)
obj_sq1 <- rep(NA, N)
eval_sq1 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], method = "squarem",
                control.method = list(K=1, version=1), control.run = list(tol=tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_sq1[i] <- end.time - start.time
    obj_sq1[i] <- fp$value.objfn
    eval_sq1[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_sq1))))
print(round(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq1[!is.na(eval_sq1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq1[!is.na(obj_sq1)], probs = c(.5, .25, .75)), 4))

##########################################
### SqS2
#############################################

time_sq2 <- rep(NA, N)
obj_sq2 <- rep(NA, N)
eval_sq2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], method = "squarem",
                control.method = list(K=1, version=2), control.run = list(tol=tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_sq2[i] <- end.time - start.time
    obj_sq2[i] <- fp$value.objfn
    eval_sq2[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_sq2))))
print(round(quantile(time_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq2[!is.na(eval_sq2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq2[!is.na(obj_sq2)], probs = c(.5, .25, .75)), 4))

##########################################
### SqS3
#############################################

time_sq3 <- rep(NA, N)
obj_sq3 <- rep(NA, N)
eval_sq3 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], method = "squarem",
                control.method = list(K=1, version=3), control.run = list(tol=tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_sq3[i] <- end.time - start.time
    obj_sq3[i] <- fp$value.objfn
    eval_sq3[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_sq3))))
print(round(quantile(time_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq3[!is.na(eval_sq3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq3[!is.na(obj_sq3)], probs = c(.5, .25, .75)), 4))

##########################################
## ZAL, q=1
##########################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], pconstr = param_constraint,
                method = "qn", control.method = list(qn=1), control.run = list(tol = tol, maxiter = 5e4))
  end.time <- Sys.time()
  if(fp$convergence){
    time_zal[i] <- end.time - start.time
    obj_zal[i] <- fp$value.objfn
    eval_zal[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_zal))))
print(round(quantile(time_zal[!is.na(time_zal)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal[!is.na(eval_zal)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal[!is.na(obj_zal)], probs = c(.5, .25, .75)), 4))

##########################################
## ZAL, q=2
##########################################

time_zal2 <- rep(NA, N)
obj_zal2 <- rep(NA, N)
eval_zal2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], pconstr = param_constraint,
                method = "qn", control.method = list(qn=2), control.run = list(tol = tol, maxiter = 5e4))
  end.time <- Sys.time()
  if(fp$convergence){
    time_zal2[i] <- end.time - start.time
    obj_zal2[i] <- fp$value.objfn
    eval_zal2[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_zal2))))
print(round(quantile(time_zal2[!is.na(time_zal2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal2[!is.na(eval_zal2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal2[!is.na(obj_zal2)], probs = c(.5, .25, .75)), 4))

##########################################
## ZAL, q=min(p, 10)
##########################################

time_zal3 <- rep(NA, N)
obj_zal3 <- rep(NA, N)
eval_zal3 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data[i,,], pconstr = param_constraint,
                method = "qn", control.method = list(qn=min(P,10)), control.run = list(tol = tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_zal3[i] <- end.time - start.time
    obj_zal3[i] <- fp$value.objfn
    eval_zal3[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_zal3))))
print(round(quantile(time_zal3[!is.na(time_zal3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal3[!is.na(eval_zal3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal3[!is.na(obj_zal3)], probs = c(.5, .25, .75)), 4))

##########################################
## DAAREM
##########################################

time_dar <- rep(NA, N)
obj_dar <- rep(NA, N)
eval_dar <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- daarem(par = start, fixptfn = update, objfn = neg.objective, n=n, dim=dim, data=data[i,,], control = list(tol = tol, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_dar[i] <- end.time - start.time
    obj_dar[i] <- -fp$value.objfn
    eval_dar[i] <- fp$fpeval
  }
}

print(paste("Number of failures:", sum(is.na(time_dar))))
print(round(quantile(time_dar, probs = c(.5, .25, .75)), 3))
print(quantile(eval_dar, probs = c(.5, .25, .75)))
print(round(quantile(obj_dar, probs = c(.5, .25, .75)), 4))


save(time_mm, time_pxem, time_bqn1, time_bqn2, time_bqn3, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, time_zal2, time_zal3, time_dar,
     eval_mm, eval_pxem, eval_bqn1, eval_bqn2, eval_bqn3, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal, eval_zal2, eval_zal3, eval_dar,
     obj_mm, obj_pxem, obj_bqn1, obj_bqn2, obj_bqn3, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, obj_zal2, obj_zal3, obj_dar, file = "Out/multiT-objects.Rdata")

load(file = "Out/multiT-objects.Rdata")

