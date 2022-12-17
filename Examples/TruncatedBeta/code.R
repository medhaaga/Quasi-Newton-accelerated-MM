###############################################
######### Truncated beta binomial #############
###############################################

rm(list = ls())
set.seed(1)
library(SQUAREM)
library(quasiNewtonMM)
library(daarem)
library(turboEM)
source("qnamm.r")
source("functions.R")
library(RColorBrewer)
p <- 2
batch <- 4
data <- matrix(c(15, 5, 2, 2, 12, 6, 7, 6, 10, 9, 2, 7, 26, 15, 3, 9), ncol = 4, nrow = 4, byrow = TRUE)
tol <- 1e-7
for (h in 1:4)
{
  freq1 <- data[h,1]
  freq2 <- data[h,2]
  freq3 <- data[h,3]
  freq4 <- data[h,4]
  freq <- freq1 + freq2 + freq3 + freq4
  start <- c(.5, 1)

  x <- seq(0.000001, 0.6, .01)
  y <- seq(0, 2, 0.01)
  z <- outer(X=x, Y=y, f, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)


  ###########################################
  ## Unaccelerated MM Algorithm
  ###########################################

  now <- start
  new <- start
  diff <- 100
  iter_mm <- 1
  chain_mm <- matrix(0, nrow = 1e7, ncol = 2)
  chain_mm[iter_mm,] <- new
  l.now <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)
  l.new <- l.now
  start.time <- Sys.time()
  while((diff > tol))
  {
    iter_mm <- iter_mm + 1
    if(iter_mm %% 1000 == 0) print(iter_mm)
    new <- update(now, batch, freq1, freq2, freq3, freq4)
    l.new <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)
    chain_mm[iter_mm,] <- new
    diff <- abs(l.new - l.now)
    now <- new
    l.now <- l.new
  }
  end.time <- Sys.time()
  time_mm <- end.time - start.time
  obj_mm <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)

  print(paste("F evals: ", iter_mm, "Time: ", round(time_mm, 3), "Negative log likelihood: ", round(obj_mm, 5)))

  pdf(file = paste("Out/beta-contour", h, "_MM.pdf", sep = ""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain_mm[1:iter_mm,1],chain_mm[1:iter_mm,2], col = c(rep(1,(iter_mm-1)), 2), pch = c(rep(1,(iter_mm-1)), 19), cex = c(rep(2,(iter_mm-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()

  ####################################
  #### J & J QN1
  #####################################

  broy_fun <- function(x, func, ...){
    return(func(x, ...) - x)
  }
  now <- as.matrix(start, ncol=1)
  new <- start
  G_now <- broy_fun(now, func = update, batch = 4, freq1=freq1, freq2=freq2,
                    freq3=freq3, freq4=freq4)
  G_new <- G_now
  qn1_chain <- matrix(now, ncol = 1)
  H <- -diag(p)
  itr <- 1
  diff <- 100
  start.time <- Sys.time()
  while(diff > tol){
    itr <- itr+1
    new <- now - H%*%G_now
    l_new <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)
    if(is.na(l_new)){
      print("Falling back to MM step")
      new <- G_now + now}
    qn1_chain <- cbind(qn1_chain, new)
    G_new <- broy_fun(new, func = update, batch = 4, freq1=freq1, freq2=freq2,
                      freq3=freq3, freq4=freq4)
    foo <- H%*%(G_new - G_now)
    H <- H + (((new-now) - foo)/as.numeric(crossprod((new-now), foo)))%*%(t(new-now)%*%H)
    diff <- norm(new-now, "2")
    now <- new
    G_now <- G_new
  }
  end.time <- Sys.time()

  time_qn1 <- end.time - start.time
  obj_qn1 <- log.likelihood(new, batch, freq1, freq2, freq3, freq4)
  iter_qn1 <- itr

  print(paste("F evals: ", iter_qn1, "Time: ", round(time_qn1, 3), "Negative log likelihood: ", round(obj_qn1, 5)))
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(qn1_chain[1,],qn1_chain[2,], col = c(rep(1,(itr-1)), 2), pch = c(rep(1,(itr-1)), 19), cex = c(rep(2,(itr-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))


  ########################################
  ## BQN, q=1
  ########################################

  start.time <- Sys.time()
  fp_bqn1 <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                 freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                 control = list(tol = tol, qn=1, objfn.inc=1e-5, step.max=1e6, maxiter = 5e4, intermed = TRUE, obj.stop=TRUE))
  end.time <- Sys.time()
  time_bqn1 <- end.time - start.time

  pdf(file = paste("Out/beta-contour", h, "_BQN1.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn1$p.inter[,1],fp_bqn1$p.inter[,2], col = c(rep(1,(fp_bqn1$iter-1)), 2), pch = c(rep(1,(fp_bqn1$iter-1)), 19), cex = c(rep(2,(fp_bqn1$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()

  print(paste("Fevals: ", fp_bqn1$fpevals, "Iterations: ", fp_bqn1$iter, "Time: ", round(time_bqn1, 3),
              "Negative log likelihood: ", round(fp_bqn1$value.objfn, 5)))

  ########################################
  ## BQN, q=2
  ########################################

  start.time <- Sys.time()
  fp_bqn2 <- BQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                 freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                 control = list(tol = tol, qn=2, objfn.inc = 1e-5, maxiter = 5e4, intermed = TRUE, obj.stop=TRUE))
  end.time <- Sys.time()
  time_bqn2 <- end.time - start.time

  print(paste("F evals: ", fp_bqn2$fpevals, "Iterations: ", fp_bqn2$iter, "Time: ", round(time_bqn2, 3),
              "Negative log likelihood: ", round(fp_bqn2$value.objfn, 5)))

  pdf(file = paste("Out/beta-contour", h, "_BQN2.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn2$p.inter[,1],fp_bqn2$p.inter[,2], col = c(rep(1,(fp_bqn2$iter-1)), 2), pch = c(rep(1,(fp_bqn2$iter-1)), 19), cex = c(rep(2,(fp_bqn2$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()


  #########################################
  ## L-BQN
  ########################################
  # heuristic fix for oscillating algo near minimum


  start.time <- Sys.time()
  fp_lbqn <- LBQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                  freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                  control = list(m=2, tol = tol, objfn.inc = 1e-5, maxiter = 5e4, intermed = TRUE, obj.stop=TRUE))
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
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=1, objfn.inc=1e-5), control.run = list(tol=tol, maxiter = 5e4, convtype = "objfn"))
  end.time <- Sys.time()
  time_sq1 <- end.time - start.time

  print(paste("F evals: ", fp_sq1$fpeval, "Iterations: ", fp_sq1$itr, "Time: ", round(time_sq1, 3), "Negative log likelihood: ", round(fp_sq1$value.objfn, 4)))


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
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=2, objfn.inc=1e-5), control.run = list(tol=tol, maxiter = 5e4, convtype = "objfn"))
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
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "squarem", control.method = list(K=1, version=3, objfn.inc=1e-5), control.run = list(tol=tol, maxiter = 5e4, convtype = "objfn"))
  end.time <- Sys.time()
  time_sq3 <- end.time - start.time

  print(paste("F evals: ", fp_sq3$fpeval, "Iterations: ", fp_sq3$itr, "Time: ", round(time_sq3, 3), "Negative log likelihood: ", round(fp_sq3$value.objfn, 4)))

  ###################################################################
  ########################### ZAL, q=1  #############################
  ###################################################################

  fp_zal_temp <- qnamm(x=start, fx_mm = update, qn=1, fx_obj = log.likelihood, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)

  pdf(file = paste("Out/beta-contour", h, "_ZAL.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal_temp$Xhist[1,],fp_zal_temp$Xhist[2,],
                                                              col = c(rep(1,(ncol(fp_zal_temp$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal_temp$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal_temp$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()

  start.time <- Sys.time()
  fp_zal <- turboem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "qn", control.method = list(qn=1), control.run = list(tol=tol, maxiter = 5e4, convtype = "objfn"))
  end.time <- Sys.time()
  time_zal <- end.time - start.time
  print(paste("F evals: ", fp_zal$fpeval, "Iterations: ", fp_zal$itr, "Time: ", round(time_zal, 3), "Negative log likelihood: ", round(fp_zal$value.objfn, 4)))

  ###################################################################
  ########################### ZAL, q=2  #############################
  ###################################################################

  fp_zal_temp2 <- qnamm(x=start, fx_mm = update, qn=1, fx_obj = log.likelihood, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)

  pdf(file = paste("Out/beta-contour", h, "_ZAL2.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal_temp2$Xhist[1,],fp_zal_temp2$Xhist[2,], col = c(rep(1,(ncol(fp_zal_temp2$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal_temp2$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal_temp2$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()

  start.time <- Sys.time()
  fp_zal2 <- turboem(par = start, fixptfn = update, objfn = log.likelihood, batch = 4, pconstr = param_constraint,
                    freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4, method = "qn", control.method = list(qn=min(p,10)), control.run = list(tol=tol, maxiter = 5e4, convtype = "objfn"))
  end.time <- Sys.time()
  time_zal2 <- end.time - start.time
  print(paste("F evals: ", fp_zal2$fpeval, "Iterations: ", fp_zal2$itr, "Time: ", round(time_zal, 3), "Negative log likelihood: ", round(fp_zal2$value.objfn, 4)))

  #################################################
  ##### DAAREM
  #################################################


  start.time <- Sys.time()
  fp_dar <- daarem(par = start, fixptfn = update, objfn = neg.objective, batch = 4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
               control = list(tol = tol, maxiter = 1e5, intermed=TRUE, mon.tol = 1e-5, convtype = "objfn"))
  end.time <- Sys.time()
  time_dar = end.time - start.time
  fp_dar$value.objfn = -fp_dar$value.objfn

  pdf(file = paste("Out/beta-contour", h, "_DAR.pdf", sep=""), height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_dar$p.intermed[,2],fp_dar$p.intermed[,3], col = c(rep(1,(nrow(fp_dar$p.intermed)-1)), 2), pch = c(rep(1,(ncol(fp_dar$p.intermed)-1)), 19), cex = c(rep(2,(ncol(fp_dar$p.intermed)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
  dev.off()

  print(paste("F evals: ", fp_dar$fpeval, "Iterations: ", fp_dar$itr, "Time: ", round(time_dar, 3), "Negative log likelihood: ", round(fp_dar$value.objfn, 4)))



  ################ Save all objects for reproducibility #################

  save(time_mm, time_bqn1, time_bqn2, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, time_zal2, time_dar, iter_mm, chain_mm,
       obj_mm, fp_bqn1, fp_bqn2, fp_lbqn, fp_sq1, fp_sq1_img, fp_sq2, fp_sq2_img, fp_sq3, fp_sq3_img, fp_zal, fp_zal2, fp_zal_temp, fp_zal_temp2, fp_dar, file = paste("Out/beta-objects", h, ".Rdata", sep=""))

}

for (h in 1:4){
  load(file = paste("Out/beta-objects", h, ".Rdata", sep=""))
  print(paste("Case = ", h))
  print(paste("MM Algo ------ F evals: ", iter_mm, "Time: ", round(time_mm, 3), "Negative log likelihood: ", round(obj_mm, 4)))
  print(paste("BQN1 --------- F evals: ", fp_bqn1$fpevals, "Iterations",  fp_bqn1$iter, "Time: ", round(time_bqn1, 3), "Negative log likelihood: ", round(fp_bqn1$value.objfn, 4)))
  print(paste("BQN2 --------- F evals: ", fp_bqn2$fpevals, "Iterations",  fp_bqn2$iter, "Time: ", round(time_bqn2, 3), "Negative log likelihood: ", round(fp_bqn2$value.objfn, 4)))
  print(paste("LBQN --------- F evals: ", fp_lbqn$fpevals, "Iterations",  fp_lbqn$iter,"Time: ", round(time_lbqn, 3), "Negative log likelihood: ", round(fp_lbqn$value.objfn, 4)))
  print(paste("SQUAREM1 ----- F evals: ", fp_sq1$fpeval, "Iterations",  fp_sq1$itr, "Time: ", round(time_sq1, 3), "Negative log likelihood: ", round(fp_sq1$value.objfn, 4)))
  print(paste("SQUAREM2 ----- F evals: ", fp_sq2$fpeval, "Iterations",  fp_sq2$itr, "Time: ", round(time_sq2, 3), "Negative log likelihood: ", round(fp_sq2$value.objfn, 4)))
  print(paste("SQUAREM3 ----- F evals: ", fp_sq3$fpeval, "Iterations",  fp_sq3$itr, "Time: ", round(time_sq3, 3), "Negative log likelihood: ", round(fp_sq3$value.objfn, 4)))
  print(paste("ZAL1 ---------- F evals: ", fp_zal$fpeval, "Iterations",  fp_zal$itr, "Time: ", round(time_zal, 3), "Negative log likelihood: ", round(fp_zal$value.objfn, 4)))
  print(paste("ZAL2 ---------- F evals: ", fp_zal2$fpeval, "Iterations",  fp_zal2$itr, "Time: ", round(time_zal2, 3), "Negative log likelihood: ", round(fp_zal2$value.objfn, 4)))
  print(paste("DAAREM ------- F evals: ", fp_dar$fpeval, "Iterations",  dim(fp_dar$p.intermed)[1], "Time: ", round(time_dar, 3), "Negative log likelihood: ", round(fp_dar$value.objfn, 4)))
}


############# Analyzing L-BQN for cases b and c #####################

h <- 2
freq1 <- data[h,1]
freq2 <- data[h,2]
freq3 <- data[h,3]
freq4 <- data[h,4]
freq <- freq1 + freq2 + freq3 + freq4
start <- c(.5, 1)
start.time <- Sys.time()
fp_lbqn1 <- LBQN(par = start, fixptfn = update, objfn = log.likelihood, batch = 4,
                freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4,
                control = list(m=2, tol = tol, objfn.inc = 1e-5, maxiter = 5e4, intermed = TRUE))
end.time <- Sys.time()
n1 <- nrow(fp_lbqn1$p.inter)
pdf(file = "Out/LBQN_case2.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(seq((n1-5e2), n1), fp_lbqn1$p.inter[(n1-5e2):n1,1], type = "l", ylab = "X", xlab = "Time")
plot(seq((n1-5e2), n1), fp_lbqn1$p.inter[(n1-5e2):n1,2], type = "l", ylab = "Y", xlab = "Time")
dev.off()
