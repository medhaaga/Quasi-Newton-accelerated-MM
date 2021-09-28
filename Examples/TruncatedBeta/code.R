###############################################
######### Truncated beta binomial #############
###############################################

  
set.seed(1)
rm(list = ls())
library(SQUAREM)
library(quasiNewtonMM)
source("qnamm.r")
source("functions.R")
library(RColorBrewer)

batch <- 4
data <- matrix(c(15, 5, 2, 2, 12, 6, 7, 6, 10, 9, 2, 7, 26, 15, 3, 9), ncol = 4, nrow = 4, byrow = TRUE)

for (h in 1:4)
{
  freq1 <- data[h,1]
  freq2 <- data[h,2]
  freq3 <- data[h,3]
  freq4 <- data[h,4]
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
  iter_mm <- 0
  chain_mm <- matrix(0, nrow = 1e7, ncol = 2)
  
  start.time <- Sys.time()
  while((diff > tol))
  {
    iter_mm <- iter_mm + 1
    if(iter_mm %% 1000 == 0) print(iter_mm)
    new <- update(now, batch, freq1, freq2, freq3, freq4)
    chain_mm[iter_mm,] <- new
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  time_mm <- end.time - start.time
  obj_mm <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)
  
  print(paste("F evals: ", iter_mm, "Time: ", round(time_mm, 3), "Negative log likelihood: ", round(obj_mm, 5)))
  
  pdf(file = paste("Out/beta-contour", h, "_MM.pdf", sep = ""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain_mm[1:iter_mm,1],chain_mm[1:iter_mm,2], col = c(rep(1,(iter_mm-1)), 2), pch = c(rep(1,(iter_mm-1)), 19), cex = c(rep(2,(iter_mm-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  
  
  ########################################
  ## BQN, q=1
  ########################################
  
  start.time <- Sys.time()
  fp_bqn1 <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                 freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                 control = list(tol = tol, qn=1, step.max=1e6, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  time_bqn1 <- end.time - start.time
  
  pdf(file = paste("Out/beta-contour", h, "_BQN1.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn1$p.inter[,1],fp_bqn1$p.inter[,2], col = c(rep(1,(fp_bqn1$iter-1)), 2), pch = c(rep(1,(fp_bqn1$iter-1)), 19), cex = c(rep(2,(fp_bqn1$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  print(paste("Fevals: ", fp_bqn1$fpevals, "Ierations: ", fp_bqn1$iter, "Time: ", round(time_bqn1, 3), "Negative log likelihood: ", round(fp_bqn1$value.objfn, 5)))
  
  ########################################
  ## BQN, q=2
  ########################################
  
  start.time <- Sys.time()
  fp_bqn2 <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                 freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                 control = list(tol = tol, qn=2, step.max=1e4, objfn.inc = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  time_bqn2 <- end.time - start.time
  
  print(paste("F evals: ", fp_bqn2$fpevals, "Iterations: ", fp_bqn2$iter, "Time: ", round(time_bqn2, 3), "Negative log likelihood: ", round(fp_bqn2$value.objfn, 5)))
  
  pdf(file = paste("Out/beta-contour", h, "_BQN2.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn2$p.inter[,1],fp_bqn2$p.inter[,2], col = c(rep(1,(fp_bqn2$iter-1)), 2), pch = c(rep(1,(fp_bqn2$iter-1)), 19), cex = c(rep(2,(fp_bqn2$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  
  #########################################
  ## L-BQN
  ########################################
  
  start.time <- Sys.time()
  fp_lbqn <- LBQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                  freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                  control = list(m=10, tol = tol, objfn.inc = .001, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  time_lbqn <- end.time - start.time
  
  print(paste("F evals: ", fp_lbqn$fpevals, "Iterations: ", fp_lbqn$iter, "Time: ", round(time_lbqn, 3), "Negative log likelihood: ", round(fp_lbqn$value.objfn, 5)))
  
  pdf(file = paste("Out/beta-contour", h, "_LBQN.pdf", sep = ""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_lbqn$p.inter[,1],fp_lbqn$p.inter[,2], col = c(rep(1,(fp_lbqn$iter-1)), 2), pch = c(rep(1,(fp_lbqn$iter-1)), 19), cex = c(rep(2,(fp_lbqn$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  
  ##########################################
  ### SqS1
  #############################################
  
  fp_sq1_img <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                        freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  
  pdf(file = paste("Out/beta-contour", h, "_SqS1.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); 
    points(fp_sq1_img$p.inter[,1],fp_sq1_img$p.inter[,2], col = c(rep(1,(fp_sq1_img$iter-1)), 2), pch = c(rep(1,(fp_sq1_img$iter-1)), 19), cex = c(rep(2,(fp_sq1_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()
  
  start.time <- Sys.time()
  fp_sq1 <- turboem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=1), control.run = list(tol=tol, maxiter = 5e4))
  end.time <- Sys.time()
  time_sq1 <- end.time - start.time
  
  print(paste("F evals: ", fp_sq1$fp_sq1eval, "Iterations: ", fp_sq1$itr, "Time: ", round(time_sq1, 3), "Negative log likelihood: ", round(fp_sq1$value.objfn, 4)))
  
  
  ##########################################
  ### SqS2
  #############################################
  
  fp_sq2_img <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                        freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  
  pdf(file = paste("Out/beta-contour", h, "_SqS2.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq2_img$p.inter[,1],fp_sq2_img$p.inter[,2], col = c(rep(1,(fp_sq2_img$iter-1)), 2), pch = c(rep(1,(fp_sq2_img$iter-1)), 19), cex = c(rep(2,(fp_sq2_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  start.time <- Sys.time()
  fp_sq2 <- turboem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=2), control.run = list(tol=tol, maxiter = 5e4))
  end.time <- Sys.time()
  time_sq2 <- end.time - start.time
  
  print(paste("F evals: ", fp_sq2$fpeval, "Iterations: ", fp_sq2$itr, "Time: ", round(time_sq2, 3), "Negative log likelihood: ", round(fp_sq2$value.objfn, 4)))
  
  
  ##########################################
  ### SqS3
  #############################################
  
  fp_sq3_img <- squarem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                        freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  
  pdf(file = paste("Out/beta-contour", h, "_SqS3.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq3_img$p.inter[,1],fp_sq3_img$p.inter[,2], col = c(rep(1,(fp_sq3_img$iter-1)), 2), pch = c(rep(1,(fp_sq3_img$iter-1)), 19), cex = c(rep(2,(fp_sq3_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  start.time <- Sys.time()
  fp_sq3 <- turboem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=3), control.run = list(tol=tol, maxiter = 5e4))
  end.time <- Sys.time()
  time_sq3 <- end.time - start.time
  
  print(paste("F evals: ", fp_sq3$fpeval, "Iterations: ", fp_sq3$itr, "Time: ", round(time_sq3, 3), "Negative log likelihood: ", round(fp_sq3$value.objfn, 4)))
  
  
  ###################################################################
  ### ZAL (the naive off-the-shelf implementation of this ###########
  ######## algorithm using turboem function fails to converge) ######
  ###################################################################
  
  fp_zal <- qnamm(x=start, fx_mm = update, qn=1, fx_obj = log.likelihood, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)
  
  pdf(file = paste("Out/beta-contour", h, "_ZAL.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal$Xhist[1,],fp_zal$Xhist[2,], col = c(rep(1,(ncol(fp_zal$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()
  
  start.time <- Sys.time()
  fp_zal <- turboem(par=start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "qn", control.method = list(qn=2))
  end.time <- Sys.time()
  time_zal <- end.time - start.time
  
  
  print(paste("F evals: ", fp_zal$fpeval, "Iterations: ", fp_zal$itr, "Time: ", round(time_zal, 3), "Negative log likelihood: ", round(fp_zal$value.objfn, 4)))
  
  
  ################ Save all objects for reproducibility ################# 
  
  save(time_mm, time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, iter_mm, chain_mm, 
       obj_mm, fp_bqn1, fp_bqn2, fp_lbqn, fp_sq1, fp_sq1_img, fp_sq2, fp_sq2_img, fp_sq3, fp_sq3_img, fp_zal, file = paste("Out/beta-objects", h, ".Rdata", sep=""))
  
  
}


