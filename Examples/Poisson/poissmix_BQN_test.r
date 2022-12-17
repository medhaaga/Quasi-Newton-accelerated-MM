rm(list = ls())
set.seed(1)
library(turboEM)
library(SQUAREM)
library(daarem)
library(quasiNewtonMM)

# EM algorithm for Poisson mixture update
poissmix_update = function(par,y) {
  # The fixed point mapping giving a single E and M step of the EM algorithm
  pnew <- rep(NA,3)
  i <- 0:(length(y)-1)
  zi <- par[1]*exp(-par[2])*(par[2]^i) / (par[1]*exp(-par[2])*(par[2]^i) + (1 - par[1])*exp(-par[3])*par[3]^i)
  pnew[1] <- sum(y*zi)/sum(y)
  pnew[2] <- sum(y*i*zi)/sum(y*zi)
  pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
  par <- pnew
  return(pnew)
}

# Negative log likelihood for Poisson mixture
poissmix_loss <- function(par,y) {
  # Objective function whose local minimum is a fixed point
  # negative log-likelihood of binary poisson mixture
  i <- 0:(length(y)-1)
  loglik <- y*log(par[1]*exp(-par[2])*(par[2]^i)/exp(lgamma(i+1)) +
                    (1 - par[1])*exp(-par[3])*(par[3]^i)/exp(lgamma(i+1)))
  return ( -sum(loglik) )
}

mpoissmix_loss <- function(par,y) {
  # Objective function whose local minimum is a fixed point
  # negative log-likelihood of binary poisson mixture
  i <- 0:(length(y)-1)
  loglik <- y*log(par[1]*exp(-par[2])*(par[2]^i)/exp(lgamma(i+1)) +
                    (1 - par[1])*exp(-par[3])*(par[3]^i)/exp(lgamma(i+1)))
  return ( sum(loglik) )
}

poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

N <- 10
p0 <- cbind(runif(N), matrix(runif(2*N, 0, 5), ncol = 2))
tol <- 1e-6
dim <- 3

#############################################
########## MM Algorithms ####################
#############################################

time_mm <- rep(0, N)
eval_mm <- rep(0,N)
obj_mm <- rep(0,N)

for (i in 1:N){
  diff <- 10
  old <- p0[i,]
  new <- old
  iter <- 0
  start <- Sys.time()
  while(diff > tol){
    iter <- iter + 1
    new <- poissmix_update(old, y=poissmix.dat$freq)
    diff <- norm(new - old, type = "2")
    old <- new
  }
  end = Sys.time()
  time_mm[i] <- end-start
  eval_mm[i] <- iter
  obj_mm[i] <- poissmix_loss(new, y=poissmix.dat$freq)
}

print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

#############################################
############## BQN, q=1  ####################
#############################################

time_bqn1 <- rep(NA, N)
obj_bqn1 <- rep(NA, N)
eval_bqn1 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.bqn1 <- BQN(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
                 control = list(qn=1, step.min=5, tol = tol, maxiter = 5e4, objfn.inc = .01))
  end = Sys.time()
  time_bqn1[i] <- end-start
  eval_bqn1[i] <- fp.bqn1$fpevals
  obj_bqn1[i] <- fp.bqn1$value.objfn
}

print(round(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)), 4))

#############################################
############## BQN, q=2  ####################
#############################################

time_bqn2 <- rep(NA, N)
obj_bqn2 <- rep(NA, N)
eval_bqn2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.bqn2 <- BQN(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
              control = list(qn=2, step.min=5, tol = tol, maxiter = 5e4, objfn.inc = .01))
  end = Sys.time()
  time_bqn2[i] <- end-start
  eval_bqn2[i] <- fp.bqn2$fpevals
  obj_bqn2[i] <- fp.bqn2$value.objfn
}

print(round(quantile(time_bqn2, probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn2, probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn2, probs = c(.5, .25, .75)), 4))

############################################
################## L-BQN ####################
#############################################

time_lbqn <- rep(NA, N)
obj_lbqn <- rep(NA, N)
eval_lbqn <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.lbqn <- LBQN(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
                  control = list(m=min(dim, 10), tol = tol, maxiter = 5e4))
  end = Sys.time()
  time_lbqn[i] <- end - start
  obj_lbqn[i] <- fp.lbqn$fpevals
  eval_lbqn[i] <- fp.lbqn$value.objfn
}

print(round(quantile(time_lbqn, probs = c(.5, .25, .75)), 3))
print(quantile(eval_lbqn, probs = c(.5, .25, .75)))
print(round(quantile(obj_lbqn, probs = c(.5, .25, .75)), 4))

#############################################
############## SQUAREM-1 ####################
#############################################

time_sq1 <- rep(NA, N)
obj_sq1 <- rep(NA, N)
eval_sq1 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.sq1 <- turboem(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                    control.method = list(K=1, version=1, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
  end = Sys.time()
  time_sq1[i] <- end - start
  eval_sq1[i] <- fp.sq1$fpeval
  obj_sq1[i] <- fp.sq1$value.objfn
}

print(round(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq1[!is.na(eval_sq1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq1[!is.na(obj_sq1)], probs = c(.5, .25, .75)), 4))

#############################################
############## SQUAREM-2 ####################
#############################################

time_sq2 <- rep(NA, N)
obj_sq2 <- rep(NA, N)
eval_sq2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.sq2 <- turboem(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                    control.method = list(K=1, version=2, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
  end = Sys.time()
  time_sq2[i] <- end - start
  eval_sq2[i] <- fp.sq2$fpeval
  obj_sq2[i] <- fp.sq2$value.objfn
}

print(round(quantile(time_sq2, probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq2, probs = c(.5, .25, .75)))
print(round(quantile(obj_sq2, probs = c(.5, .25, .75)), 4))

#############################################
############## SQUAREM-3 ####################
#############################################

time_sq3 <- rep(NA, N)
obj_sq3 <- rep(NA, N)
eval_sq3 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.sq3 <- turboem(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                    control.method = list(K=1, version=3, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
  end = Sys.time()
  time_sq3[i] <- end - start
  obj_sq3[i] <- fp.sq2$fpeval
  eval_sq3[i] <- fp.sq2$value.objfn
}

print(round(quantile(time_sq3, probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq3, probs = c(.5, .25, .75)))
print(round(quantile(obj_sq3, probs = c(.5, .25, .75)), 4))

#############################################
############## ZAL, q=1 ####################
#############################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.zal1 <- turboem(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "qn",
                    control.method = list(qn=1), control.run = list(tol=tol, maxiter = 5e4))
  end = Sys.time()
  time_zal[i] <- end - start
  eval_zal[i] <- fp.zal1$fpeval
  obj_zal[i] <- fp.zal1$value.objfn
}

sum(is.na(eval_zal))
print(round(quantile(time_zal[!is.na(time_zal)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal[!is.na(eval_zal)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal[!is.na(obj_zal)], probs = c(.5, .25, .75)), 4))

#############################################
############## ZAL, q=2 ####################
#############################################

time_zal2 <- rep(NA, N)
obj_zal2 <- rep(NA, N)
eval_zal2 <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.zal2 <- turboem(par = p0[i,], fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "qn",
                     control.method = list(qn=1), control.run = list(tol=tol, maxiter = 5e4))
  end = Sys.time()
  time_zal2[i] <- end - start
  eval_zal2[i] <- fp.zal2$fpeval
  obj_zal2[i] <- fp.zal2$value.objfn
}

sum(is.na(eval_zal2))
print(round(quantile(time_zal2[!is.na(time_zal2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal2[!is.na(eval_zal2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal2[!is.na(obj_zal2)], probs = c(.5, .25, .75)), 4))

#############################################
################# DAAREM ####################
#############################################

time_dar <- rep(NA, N)
obj_dar <- rep(NA, N)
eval_dar <- rep(NA, N)

for (i in 1:N){
  print(i)
  start <- Sys.time()
  fp.da <- daarem(par = p0[i,], fixptfn = poissmix_update, objfn = mpoissmix_loss, y=poissmix.dat$freq, control = list(tol = tol, maxiter = 5e4, mon.tol = 1e-5, intermed = TRUE))
  end = Sys.time()
  time_dar[i] <- end - start
  obj_dar[i] <- fp.da$fpeval
  eval_dar[i] <- fp.da$value.objfn
}

print(round(quantile(time_dar, probs = c(.5, .25, .75)), 3))
print(quantile(eval_dar, probs = c(.5, .25, .75)))
print(round(quantile(obj_dar, probs = c(.5, .25, .75)), 4))

save(time_mm, time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, time_zal2, time_dar,
     eval_mm, eval_bqn1, eval_bqn2, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal, eval_zal2, eval_dar,
     obj_mm, obj_bqn1, obj_bqn2, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, obj_zal2, obj_dar, file = "Out/poisson-objects.Rdata")



