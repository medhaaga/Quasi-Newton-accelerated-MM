rm(list = ls())
set.seed(10)
library(turboEM)
library(SQUAREM)
library(daarem)
library(quasiNewtonMM)

# EM algorithm for Poisson mixture update
poissmix_update = function(par,y) {
  # The fixed point mapping giving a single E and M step of the EM algorithm
  pnew <- rep(NA,3)
  i <- 0:(length(y)-1)
  zi <- par[1]*exp(-par[2])*par[2]^i / (par[1]*exp(-par[2])*par[2]^i + (1 - par[1])*exp(-par[3])*par[3]^i)
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
  loglik <- y*log(par[1]*exp(-par[2])*par[2]^i/exp(lgamma(i+1)) +
                    (1 - par[1])*exp(-par[3])*par[3]^i/exp(lgamma(i+1)))
  return ( -sum(loglik) )
}

mpoissmix_loss <- function(par,y) {
  # Objective function whose local minimum is a fixed point
  # negative log-likelihood of binary poisson mixture
  i <- 0:(length(y)-1)
  loglik <- y*log(par[1]*exp(-par[2])*par[2]^i/exp(lgamma(i+1)) +
                    (1 - par[1])*exp(-par[3])*par[3]^i/exp(lgamma(i+1)))
  return ( sum(loglik) )
}

poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

p0 <- c(runif(1), runif(2, 0, 5))
tol <- 1e-6
dim <- 3

#############################################
########## MM Algorithms ####################
#############################################

diff <- 10
old <- p0
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
print(paste("Time: ", end - start, "F evals: ", iter, "Objective: ", poissmix_loss(new, y=poissmix.dat$freq)))


#############################################
############## BQN, q=2  ####################
#############################################

start <- Sys.time()
fp.bqn1 <- BQN(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
              control = list(qn=2, step.min=1, tol = tol, maxiter = 5e4, objfn.inc = .01))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.bqn1$fpevals, "Objective: ", fp.bqn1$value.objfn))

#############################################
############## BQN, q=3  ####################
#############################################

start <- Sys.time()
fp.bqn2 <- BQN(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
              control = list(qn=2, step.min=1, tol = tol, maxiter = 5e4, objfn.inc = .01, intermed = TRUE))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.bqn2$fpevals, "Objective: ", fp.bqn2$value.objfn))

plot_ly(x=fp.bqn2$p.inter[,1], y=fp.bqn2$p.inter[,2], z=fp.bqn2$p.inter[,3])

#############################################
################## L-BQN ####################
#############################################

start <- Sys.time()
fp.lbqn <- LBQN(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq,
                control = list(m=min(dim, 10), tol = tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.lbqn$fpevals, "Objective: ", fp.lbqn$value.objfn))

#############################################
############## SQUAREM-1 ####################
#############################################

start <- Sys.time()
fp.sq1 <- turboem(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                  control.method = list(K=1, version=1, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.sq1$fpeval, "Objective: ", fp.sq1$value.objfn))

#############################################
############## SQUAREM-2 ####################
#############################################

start <- Sys.time()
fp.sq2 <- turboem(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                 control.method = list(K=2, version=1, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.sq2$fpeval, "Objective: ", fp.sq2$value.objfn))

#############################################
############## SQUAREM-3 ####################
#############################################

start <- Sys.time()
fp.sq3 <- turboem(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "squarem",
                  control.method = list(K=3, version=1, objfn.inc=1e-2), control.run = list(tol=tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.sq3$fpeval, "Objective: ", fp.sq3$value.objfn))

#############################################
############## ZAL, q=2 ####################
#############################################

start <- Sys.time()
fp.zal1 <- turboem(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "qn",
                  control.method = list(qn=2), control.run = list(tol=tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.zal1$fpeval, "Objective: ", fp.zal1$value.objfn))

#############################################
############## ZAL, q=3 ####################
#############################################

start <- Sys.time()
fp.zal2 <- turboem(par = p0, fixptfn = poissmix_update, objfn = poissmix_loss, y=poissmix.dat$freq, method = "qn",
                   control.method = list(qn=1), control.run = list(tol=tol, maxiter = 5e4))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.zal2$fpeval, "Objective: ", fp.zal2$value.objfn))

#############################################
################# DAAREM ####################
#############################################

start <- Sys.time()
fp.da <- daarem(par = p0, fixptfn = poissmix_update, objfn = mpoissmix_loss, y=poissmix.dat$freq, control = list(tol = tol, maxiter = 5e4, mon.tol = 1e-5, intermed = TRUE))
end = Sys.time()
print(paste("Time: ", end - start, "F evals: ", fp.da$fpeval, "Objective: ", fp.da$value.objfn))

