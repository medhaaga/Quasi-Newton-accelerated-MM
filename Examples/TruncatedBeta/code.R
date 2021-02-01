set.seed(1)
library(pracma)
library(SQUAREM)
library(quasiNewtonMM)
library(MASS)
source("qnamm.r")
library(RColorBrewer)
library(graphics)


log.likelihood <- function(x, batch, freq1, freq2, freq3, freq4){

  pi <- x[1]
  alpha <- x[2]
  freq <- freq1 + freq2 + freq3 + freq4

  distri <- function(x, batch = 4, pi, alpha){
    prod1 <- 1
    prod2 <- 1
    prod3 <- 1
    if (x != 0){
      for (j in 0:(x-1)){
        prod1 <- prod1*(pi + j*alpha)
      }
    }
    if(x != batch){
      for (k in 0:(batch-x-1)){
        prod2 <- prod2*(1 - pi + k*alpha)
      }
    }
    for (l in 0:(batch-1)){
      prod3 <- prod3*(1 + l*alpha)
    }
    return ((choose(batch, x) * prod1 * prod2)/prod3)
  }


  denominator <- freq*log(1 - distri(x = 0, batch = 4, pi, alpha))
  dens1 <- freq1*log(distri(x=1, batch=4, pi=pi, alpha))
  dens2 <- freq2*log(distri(x=2, batch=4, pi, alpha))
  dens3 <- freq3*log(distri(x=3, batch=4, pi, alpha))
  dens4 <- freq4*log(distri(x=4, batch=4, pi, alpha))

  return (-dens1 - dens2 - dens3 - dens4 + denominator)
}

update <- function(now, batch = 4, freq1, freq2, freq3, freq4){

  new <- now
  pi <- now[1]
  alpha <- now[2]
  freq <- freq1 + freq2 + freq3 + freq4
  prod2 <- 1
  prod3 <- 1
  for (k in 0:(batch-1)){
    prod2 <- prod2*(1 - pi + k*alpha)
    prod3 <- prod3*(1 + k*alpha)
  }
  foo <- prod2/prod3
  foo <- foo/(1 - foo)
  s1 <- c(freq, freq2 + freq3 + freq4, freq3 + freq4, freq4)
  s2 <- c(freq1 + freq2 + freq3 + freq*foo, freq1 + freq2 + freq*foo, freq1 + freq*foo, freq*foo)
  r <- rep(freq*(1+foo), 4)
  num1 <- 0
  num2 <- 0
  denom1 <- 0
  denom2 <- 0
  for (k in 0:3){
    num1 <- num1 + (((s1[k+1]*k*alpha)/(pi + k*alpha)) + ((s2[k+1]*k*alpha)/(1 - pi + k*alpha)))
    num2 <- num2 + ((s1[k+1]*pi)/(pi + k*alpha))
    denom1 <- denom1 + ((r[k+1]*k)/(1 + k*alpha))
    denom2 <- denom2 + (((s1[k+1]*pi) / (pi + k*alpha)) + ((s2[k+1]*(1 - pi)) / (1 - pi + k*alpha)))
  }
  new[2] <- num1/denom1
  new[1] <- num2/denom2
  return (new)
}


f <- function(pi, alpha, batch, freq1, freq2, freq3, freq4){

  freq <- freq1 + freq2 + freq3 + freq4

  distri <- function(x, batch = 4, pi, alpha){
    prod1 <- 1
    prod2 <- 1
    prod3 <- 1
    if (x != 0){
      for (j in 0:(x-1)){
        prod1 <- prod1*(pi + j*alpha)
      }
    }
    if(x != batch){
      for (k in 0:(batch-x-1)){
        prod2 <- prod2*(1 - pi + k*alpha)
      }
    }
    for (l in 0:(batch-1)){
      prod3 <- prod3*(1 + l*alpha)
    }
    return ((choose(batch, x) * prod1 * prod2)/prod3)
  }


  denominator <- freq*log(1 - distri(x = 0, batch = 4, pi, alpha))
  dens1 <- freq1*log(distri(x=1, batch=4, pi=pi, alpha))
  dens2 <- freq2*log(distri(x=2, batch=4, pi, alpha))
  dens3 <- freq3*log(distri(x=3, batch=4, pi, alpha))
  dens4 <- freq4*log(distri(x=4, batch=4, pi, alpha))

  return (-dens1 - dens2 - dens3 - dens4 + denominator)
}


batch <- 4
freq1 <- 15
freq2 <- 5
freq3 <- 2
freq4 <- 2
freq <- freq1 + freq2 + freq3 + freq4
start <- c(.5, 1)
tol <- 1e-7

x <- seq(0.000001, 0.6, .01)
y <- seq(0, 2, 0.01)
z <- outer(X=x, Y=y, f, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)


###########################################
## MM Algorithm
###########################################

now <- start
new <- start
diff <- 100
iter <- 0
start.time <- Sys.time()
chain <- matrix(0, nrow = 1e7, ncol = 2)
while((diff > tol))
{
  iter <- iter + 1
  if(iter %% 1000 == 0) print(iter)
  new <- update(now, batch, freq1, freq2, freq3, freq4)
  chain[iter,] <- new
  diff <- sqrt(crossprod(new-now))
  now <- new
}
end.time <- Sys.time()
pdf(file = "Out/beta-contour_MM.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); 
  points(chain[1:iter,1],chain[1:iter,2], col = c(rep(1,(iter-1)), 2), pch = c(rep(1,(iter-1)), 19), cex = c(rep(1.5,(iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(paste("Fevals: ", iter, "theta: ", new, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", log.likelihood(new, batch, freq1, freq2, freq3, freq4)))


########################################
## BQN, q=1
########################################

start.time <- Sys.time()
fp <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
           freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
           control = list(tol = tol, qn=1, step.max=1000000, objfn.inc = 1, maxiter = 1e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_BQN1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1.5,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))

########################################
## BQN, q=2
########################################

start.time <- Sys.time()
fp <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
           freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
           control = list(tol = tol, qn=2, step.max=1000000, objfn.inc = 1, maxiter = 1e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_BQN2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1.5,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))


  ########################################
## L-BQN
########################################

start.time <- Sys.time()
fp <- LBQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
           freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
           control = list(m=10, tol = tol, objfn.inc = .001, maxiter = 1e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_LBFGS.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))


##########################################
### SqS1
#############################################

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
              freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_SqS1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); 
  points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1.5,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))


##########################################
### SqS2
#############################################

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
              freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_SqS2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1.5,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))


##########################################
### SqS3
#############################################

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
              freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
end.time <- Sys.time()

pdf(file = "Out/beta-contour_SqS3.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2], col = c(rep(1,(fp$iter-1)), 2), pch = c(rep(1,(fp$iter-1)), 19), cex = c(rep(1,(fp$iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(paste("Fevals: ", fp$fpevals, "Ierations: ", fp$iter, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$value.objfn))


######################################
### ZAL
########################################

start.time <- Sys.time()
fp <- qnamm(x=start, fx_mm=update, qn=1, fx_obj=log.likelihood, max_iter=5e4, tol=tol, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)
end.time <- Sys.time()

iter <- dim(fp$Xhist)[2]
pdf(file = "Out/beta-contour_ZAL1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$Xhist[1,],fp$Xhist[2,], col = c(rep(1,(iter-1)), 2), pch = c(rep(1,(iter-1)), 19), cex = c(rep(1.5,(iter-1)), 2))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(paste("Fevals: ", fp$fevals, "Ierations: ", fp$accept + fp$reject, "Time: ", end.time-start.time, 
            "Negative log likelihood: ", fp$objective))

