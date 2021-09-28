#################################################
######## Generalised Eigenvalues ################
#################################################
rm(list = ls())
library(pracma)
library(turboEM)
library(quasiNewtonMM)
source("functions.R")

set.seed(1)
dim <- 100
C <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
D <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
A <- C + t(C)
B <- D %*% t(D)

N <- 10
start_rep <- matrix(rnorm(N*dim, mean = 0, sd = 100), nrow = N, ncol = dim)
tol <- 1e-7
D <- c("descent", "ascent")

for (d in 1:2)
{
  dir <- D[d]
  ###########################################
  ## Unaccelerated MM Algorithm
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
      new <- update(now, A, B, dir = dir)
      diff <- sqrt(crossprod(new-now))
      now <- new
    }
    end.time <- Sys.time()
    
    time_mm[i] <- end.time - start.time
    obj_mm[i] <- rayleigh(new, A, B, dir)
    eval_mm[i] <- iter
  }
  
  print(quantile(time_mm, probs = c(.5, .25, .75)))
  print(quantile(eval_mm, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_mm, probs = c(.5, .25, .75)))
  
  
  ########################################
  ## BQN, q=1
  ########################################
  
  time_bqn1 <- rep(NA, N)
  obj_bqn1 <- rep(NA, N)
  eval_bqn1 <- rep(NA, N)
  
  for (i in 1:N){
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(qn=1, tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_bqn1[i] <- end.time - start.time
      obj_bqn1[i] <- fp$value.objfn
      eval_bqn1[i] <- fp$fpevals
    }
  }
  
  print(paste("Number of failures: ", sum(is.na(time_bqn1))))
  print(quantile(time_bqn1, probs = c(.5, .25, .75)))
  print(quantile(eval_bqn1, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_bqn1, probs = c(.5, .25, .75)))
  
  #######################################
  ## L-BQN
  ########################################
  
  time_lbqn <- rep(NA, N)
  obj_lbqn <- rep(NA, N)
  eval_lbqn <- rep(NA, N)
  
  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- LBQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(m=10, tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_lbqn[i] <- end.time - start.time
      obj_lbqn[i] <- fp$value.objfn
      eval_lbqn[i] <- fp$fpevals
    }
    
  }
  print(paste("Number of failures: ", sum(is.na(time_lbqn))))
  print(quantile(time_lbqn, probs = c(.5, .25, .75)))
  print(quantile(eval_lbqn, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_lbqn, probs = c(.5, .25, .75)))
  
  ##########################################
  ### SQUAREM v1
  #############################################
  
  time_sq1 <- rep(NA, N)
  obj_sq1 <- rep(NA, N)
  eval_sq1 <- rep(NA, N)
  
  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, method = "squarem", control.method = list(K=1, version=1), control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_sq1[i] <- end.time - start.time
      obj_sq1[i] <- fp$value.objfn
      eval_sq1[i] <- fp$fpeval
    }
  }
  
  print(paste("Number of failures: ", sum(is.na(time_sq1))))
  print(quantile(time_sq1, probs = c(.5, .25, .75)))
  print(quantile(eval_sq1, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq1, probs = c(.5, .25, .75)))
  
  ##########################################
  ### SQUAREM v2
  #############################################
  
  time_sq2 <- rep(NA, N)
  obj_sq2 <- rep(NA, N)
  eval_sq2 <- rep(NA, N)
  
  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, method = "squarem", control.method = list(K=1, version=2), control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_sq2[i] <- end.time - start.time
      obj_sq2[i] <- fp$value.objfn
      eval_sq2[i] <- fp$fpeval
    }
  }
  
  print(paste("Number of failures: ", sum(is.na(time_sq2))))
  print(quantile(time_sq2, probs = c(.5, .25, .75)))
  print(quantile(eval_sq2, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq2, probs = c(.5, .25, .75)))
  
  
  ##########################################
  ### SQUAREM v3
  #############################################
  
  time_sq3 <- rep(NA, N)
  obj_sq3 <- rep(NA, N)
  eval_sq3 <- rep(NA, N)
  
  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, method = "squarem", control.method = list(K=1, version=3), control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_sq3[i] <- end.time - start.time
      obj_sq3[i] <- fp$value.objfn
      eval_sq3[i] <- fp$fpeval
    }
  }
  
  print(paste("Number of failures: ", sum(is.na(time_sq3))))
  print(quantile(time_sq3, probs = c(.5, .25, .75)))
  print(quantile(eval_sq3, probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq3, probs = c(.5, .25, .75)))
  
  ##########################################
  ## ZAL for q=2
  ##########################################
  
  time_zal <- rep(NA, N)
  obj_zal <- rep(NA, N)
  eval_zal <- rep(NA, N)
  
  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, method = "qn", control.method = list(qn=2), control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    
    if(fp$convergence){
      time_zal[i] <- end.time - start.time
      obj_zal[i] <- fp$value.objfn
      eval_zal[i] <- fp$fpeval
    }
    
  }
  
  print(paste("Number of failures: ", sum(is.na(time_zal))))
  print(quantile(time_zal))
  print(quantile(eval_zal))
  print(quantile((-1)^(d+1)*obj_zal))
  
  ##############################################
  
  save(time_mm, time_bqn1, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal,
       eval_mm, eval_bqn1, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal,
       obj_mm, obj_bqn1, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = paste("Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))

  
}

###########################################################
###### Creating scatter plots #############################
###########################################################

for ( d in 1:2)
{
  dir <- D[d]
  
  load(paste("Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))
  
  time_range <- range(time_mm, time_bqn1, time_lbqn, time_sq3, time_zal)
  eval_range <- range(eval_mm, eval_bqn1, eval_lbqn, eval_sq3,  eval_zal)
  obj_range <- (-1)^(d+1)*range(obj_mm, obj_bqn1, obj_lbqn, obj_sq3, obj_zal)
  
  pdf(file = paste("Out/eigen-objVSeval_", dir, ".pdf", sep = ""), height = 5, width = 6)
  par(mar=c(5, 4, 4, 8), xpd = TRUE)
  plot(eval_bqn1, (-1)^(d+1)*obj_bqn1, pch=19, cex=1.2, col  ="red", xlim = eval_range, ylim = obj_range, ylab = "Objective Value", xlab = "Number of F evaluations")
  points(eval_lbqn, (-1)^(d+1)*obj_lbqn, pch=19, cex=1.2, col = "blue")
  points(eval_sq3, (-1)^(d+1)*obj_sq3, pch=19, cex=1.2, col = "pink")
  points(eval_zal, (-1)^(d+1)*obj_zal, pch=19, cex=1.2, col = "orange")
  points(eval_mm, (-1)^(d+1)*obj_mm, pch=19, cex=1.2, col = "black")
  legend("topright", inset = c(-.4, 0), col = c("black", "red", "blue", "pink", "orange"), pch=19, cex=1.2, legend  =c("MM", "BQN, q=1", "L-BQN", "SQUAREM v3", "ZAL, q=2"))
  dev.off()
  
  pdf(file = paste("Out/eigen-objVStime_", dir, ".pdf", sep = ""), height = 5, width = 6)
  par(mar=c(5, 4, 4, 8), xpd = TRUE)
  plot(time_bqn1, (-1)^(d+1)*obj_bqn1, pch=19, cex=1.2, col  ="red", xlim = time_range, ylim = obj_range, ylab = "Objective Value", xlab = "Time (s)")
  points(time_lbqn, (-1)^(d+1)*obj_lbqn, pch=19, cex=1.2, col = "blue")
  points(time_sq3, (-1)^(d+1)*obj_sq3, pch=19, cex=1.2, col = "pink")
  points(time_zal, (-1)^(d+1)*obj_zal, pch=19, cex=1.2, col = "orange")
  points(time_mm, (-1)^(d+1)*obj_mm, pch=19, cex=1.2, col = "black")
  legend("topright", inset = c(-.4, 0), col = c("black", "red", "blue", "pink", "orange"), pch=19, cex=1.2, legend  =c("MM", "BQN, q=1", "L-BQN", "SQUAREM v3", "ZAL, q=2"))
  dev.off()
}

