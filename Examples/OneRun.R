set.seed(1)
rm(list = ls())
library(RColorBrewer)
source("TruncatedBeta/functions.R")


########################################
###### Example 1: Table-1 ##############
########################################

load(file = "Quadratic/Out/quad-objects_sq1e2.Rdata")

print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

print(round(quantile(time_bqn1, c(.5, .25, .75)), 3))
print(quantile(eval_bqn1, c(.5, .25, .75)))
print(round(quantile(obj_bqn1, c(.5, .25, .75)), 4))

print(round(quantile(time_bqn2, c(.5, .25, .75)), 3))
print(quantile(eval_bqn2, c(.5, .25, .75)))
print(round(quantile(obj_bqn2, c(.5, .25, .75)), 4))

print(round(quantile(time_lbqn, c(.5, .25, .75)), 3))
print(quantile(eval_lbqn, c(.5, .25, .75)))
print(round(quantile(obj_lbqn, c(.5, .25, .75)), 4))

print(round(quantile(time_sq1, c(.5, .25, .75)), 3))
print(quantile(eval_sq1, c(.5, .25, .75)))
print(round(quantile(obj_sq1, c(.5, .25, .75)), 4))

print(round(quantile(time_sq2, c(.5, .25, .75)), 3))
print(quantile(eval_sq2, c(.5, .25, .75)))
print(round(quantile(obj_sq2, c(.5, .25, .75)), 4))

print(round(quantile(time_sq3, c(.5, .25, .75)), 3))
print(quantile(eval_sq3, c(.5, .25, .75)))
print(round(quantile(obj_sq3, c(.5, .25, .75)), 4))

print(round(quantile(time_zal, c(.5, .25, .75)), 3))
print(quantile(eval_zal, c(.5, .25, .75)))
print(round(quantile(obj_zal, c(.5, .25, .75)), 4))


########################################
###### Example 1: Figure-2 #############
########################################

load(file = "Quadratic/Out/quad-objects_sq1e2.Rdata")

df1 <- data.frame( "B1" = eval_bqn1, "B2" = eval_bqn2, "L-B" = eval_lbqn, "Sq" = eval_sq3, "ZAL" = eval_zal)

df2 <- data.frame("B1" = time_bqn1, "B2" = time_bqn2, "L-B" = time_lbqn, "Sq" = time_sq3, "ZAL" = time_zal)


pdf(file = "Quadratic/Out/quad-boxplot_sd1e2.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
boxplot(df1, xlab = "Acceleration Method", ylab = "No. of evaluations")
boxplot(df2, xlab = "Acceleration Method", ylab = "Time")
dev.off()

########################################
###### Example 2: Table-2 #############
########################################

load(file = "TruncatedBeta/Out/beta-objects.Rdata")

print(paste("F evals: ", iter_mm, "Time: ", round(time_mm, 3), "Negative log likelihood: ", round(obj_mm, 5)))
print(paste("Fevals: ", fp_bqn1$fpevals, "Ierations: ", fp_bqn1$iter, "Time: ", round(time_bqn1, 3), "Negative log likelihood: ", round(fp_bqn1$value.objfn, 5)))
print(paste("F evals: ", fp_bqn2$fpevals, "Iterations: ", fp_bqn2$iter, "Time: ", round(time_bqn2, 3), "Negative log likelihood: ", round(fp_bqn2$value.objfn, 5)))
print(paste("F evals: ", fp_lbqn$fpevals, "Iterations: ", fp_lbqn$iter, "Time: ", round(time_lbqn, 3), "Negative log likelihood: ", round(fp_lbqn$value.objfn, 5)))
print(paste("F evals: ", fp_sq1$fp_sq1eval, "Iterations: ", fp_sq1$itr, "Time: ", round(time_sq1, 3), "Negative log likelihood: ", round(fp_sq1$value.objfn, 4)))
print(paste("F evals: ", fp_sq2$fpeval, "Iterations: ", fp_sq2$itr, "Time: ", round(time_sq2, 3), "Negative log likelihood: ", round(fp_sq2$value.objfn, 4)))
print(paste("F evals: ", fp_sq3$fpeval, "Iterations: ", fp_sq3$itr, "Time: ", round(time_sq3, 3), "Negative log likelihood: ", round(fp_sq3$value.objfn, 4)))
print(paste("F evals: ", fp_zal$fevals, "Iterations: ", fp_zal$accept + fp_zal$reject, "Time: ", round(time_zal, 3), "Negative log likelihood: ", round(fp_zal$objective, 4)))

########################################
###### Example 2: Figure-3 #############
########################################

load(file = "TruncatedBeta/Out/beta-objects.Rdata")

batch <- 4
freq1 <- 15
freq2 <- 5
freq3 <- 2
freq4 <- 2
x <- seq(0.000001, 0.6, .01)
y <- seq(0, 2, 0.01)
z <- outer(X=x, Y=y, f, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)


pdf(file = "TruncatedBeta/Out/beta-contour_MM.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain_mm[1:iter_mm,1],chain_mm[1:iter_mm,2], col = c(rep(1,(iter_mm-1)), 2), pch = c(rep(1,(iter_mm-1)), 19), cex = c(rep(2,(iter_mm-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_BQN1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn1$p.inter[,1],fp_bqn1$p.inter[,2], col = c(rep(1,(fp_bqn1$iter-1)), 2), pch = c(rep(1,(fp_bqn1$iter-1)), 19), cex = c(rep(2,(fp_bqn1$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_BQN2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn2$p.inter[,1],fp_bqn2$p.inter[,2], col = c(rep(1,(fp_bqn2$iter-1)), 2), pch = c(rep(1,(fp_bqn2$iter-1)), 19), cex = c(rep(2,(fp_bqn2$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_LBQN.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_lbqn$p.inter[,1],fp_lbqn$p.inter[,2], col = c(rep(1,(fp_lbqn$iter-1)), 2), pch = c(rep(1,(fp_lbqn$iter-1)), 19), cex = c(rep(2,(fp_lbqn$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_SqS1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); 
  points(fp_sq1_img$p.inter[,1],fp_sq1_img$p.inter[,2], col = c(rep(1,(fp_sq1_img$iter-1)), 2), pch = c(rep(1,(fp_sq1_img$iter-1)), 19), cex = c(rep(2,(fp_sq1_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_SqS2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq2_img$p.inter[,1],fp_sq2_img$p.inter[,2], col = c(rep(1,(fp_sq2_img$iter-1)), 2), pch = c(rep(1,(fp_sq2_img$iter-1)), 19), cex = c(rep(2,(fp_sq2_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_SqS3.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq3_img$p.inter[,1],fp_sq3_img$p.inter[,2], col = c(rep(1,(fp_sq3_img$iter-1)), 2), pch = c(rep(1,(fp_sq3_img$iter-1)), 19), cex = c(rep(2,(fp_sq3_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour_ZAL.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal$Xhist[1,],fp_zal$Xhist[2,], col = c(rep(1,(ncol(fp_zal$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()


########################################
###### Example 3: Table-3 #############
########################################

D = c("descent", "ascent")

for (d in 1:2)
{
  dir <- D[d]
  load(paste("GenEigen/Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))
  
  print(paste(round(time_mm[1], 3), eval_mm[1], round((-1)^(d+1)*obj_mm[1], 4)))
  print(paste(round(time_bqn1[1], 3), eval_bqn1[1], round((-1)^(d+1)*obj_bqn1[1], 4)))
  print(paste(round(time_lbqn[1], 3), eval_lbqn[1], round((-1)^(d+1)*obj_lbqn[1], 4)))
  print(paste(round(time_sq1[1], 3), eval_sq1[1], round((-1)^(d+1)*obj_sq1[1], 4)))
  print(paste(round(time_sq2[1], 3), eval_sq2[1], round((-1)^(d+1)*obj_sq2[1], 4)))
  print(paste(round(time_sq3[1], 3), eval_sq3[1], round((-1)^(d+1)*obj_sq3[1], 4)))
  print(paste(round(time_zal[1], 3), eval_zal[1], round((-1)^(d+1)*obj_zal[1], 4)))
}


########################################
###### Example 3: Figure - 4 #############
########################################

D = c("descent", "ascent")

for ( d in 1:2)
{
  dir <- D[d]
  
  load(paste("GenEigen/Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))
  
  time_range <- range(time_mm, time_bqn1, time_lbqn, time_sq3, time_zal)
  eval_range <- range(eval_mm, eval_bqn1, eval_lbqn, eval_sq3,  eval_zal)
  obj_range <- (-1)^(d+1)*range(obj_mm, obj_bqn1, obj_lbqn, obj_sq3, obj_zal)
  
  pdf(file = paste("GenEigen/Out/eigen-objVSeval_", dir, ".pdf", sep = ""), height = 5, width = 6)
  par(mar=c(5, 4, 4, 8), xpd = TRUE)
  plot(eval_bqn1, (-1)^(d+1)*obj_bqn1, pch=19, col  ="red", xlim = eval_range, ylim = obj_range, ylab = "Objective Value", xlab = "Number of F evaluations")
  points(eval_lbqn, (-1)^(d+1)*obj_lbqn, pch=19, col = "blue")
  points(eval_sq3, (-1)^(d+1)*obj_sq3, pch=19, col = "pink")
  points(eval_zal, (-1)^(d+1)*obj_zal, pch=19, col = "orange")
  points(eval_mm, (-1)^(d+1)*obj_mm, pch=19, col = "black")
  legend("topright", inset = c(-.4, 0), col = c("black", "red", "blue", "pink", "orange"), pch=19, legend  =c("MM", "BQN, q=1", "L-BQN", "SQUAREM v3", "ZAL, q=2"))
  dev.off()
  
  pdf(file = paste("GenEigen/Out/eigen-objVStime_", dir, ".pdf", sep = ""), height = 5, width = 6)
  par(mar=c(5, 4, 4, 8), xpd = TRUE)
  plot(time_bqn1, (-1)^(d+1)*obj_bqn1, pch=19, col  ="red", xlim = time_range, ylim = obj_range, ylab = "Objective Value", xlab = "Time (s)")
  points(time_lbqn, (-1)^(d+1)*obj_lbqn, pch=19, col = "blue")
  points(time_sq3, (-1)^(d+1)*obj_sq3, pch=19, col = "pink")
  points(time_zal, (-1)^(d+1)*obj_zal, pch=19, col = "orange")
  points(time_mm, (-1)^(d+1)*obj_mm, pch=19, col = "black")
  legend("topright", inset = c(-.4, 0), col = c("black", "red", "blue", "pink", "orange"), pch=19, legend  =c("MM", "BQN, q=1", "L-BQN", "SQUAREM v3", "ZAL, q=2"))
  dev.off()
}

########################################
###### Example 4: Table - 4 #############
########################################

load(file = "MultiT/Out/multiT-objects.Rdata")

print(paste("Number of failures:", sum(is.na(time_mm))))
print(round(time_mm, 3))
print(eval_mm)
print(round(obj_mm, 4))

print(paste("Number of failures:", sum(is.na(time_pxem))))
print(round(time_pxem[!is.na(time_pxem)], 3))
print(eval_pxem[!is.na(eval_pxem)])
print(round(obj_pxem[!is.na(obj_pxem)], 4))

print(paste("Number of failures:", sum(is.na(time_bqn1))))
print(round(time_bqn1[!is.na(time_bqn1)], 3))
print(eval_bqn1[!is.na(eval_bqn1)])
print(round(obj_bqn1[!is.na(obj_bqn1)], 4))

print(paste("Number of failures:", sum(is.na(time_bqn2))))
print(round(time_bqn2[!is.na(time_bqn2)], 3))
print(eval_bqn2[!is.na(eval_bqn2)])
print(round(obj_bqn2[!is.na(obj_bqn2)], 4))

print(paste("Number of failures:", sum(is.na(time_lbqn))))
print(round(time_lbqn[!is.na(time_lbqn)], 3))
print(eval_lbqn[!is.na(eval_lbqn)])
print(round(obj_lbqn[!is.na(obj_lbqn)], 4))

print(paste("Number of failures:", sum(is.na(time_sq1))))
print(round(time_sq1[!is.na(time_sq1)], 3))
print(eval_sq1[!is.na(eval_sq1)])
print(round(obj_sq1[!is.na(obj_sq1)], 4))

print(paste("Number of failures:", sum(is.na(time_sq2))))
print(round(time_sq2[!is.na(time_sq2)], 3))
print(eval_sq2[!is.na(eval_sq2)])
print(round(obj_sq2[!is.na(obj_sq2)], 4))

print(paste("Number of failures:", sum(is.na(time_sq3))))
print(round(time_sq3[!is.na(time_sq3)], 3))
print(eval_sq3[!is.na(eval_sq3)])
print(round(obj_sq3[!is.na(obj_sq3)], 4))

print(paste("Number of failures:", sum(is.na(time_zal))))
print(round(time_zal[!is.na(time_zal)], 3))
print(eval_zal[!is.na(eval_zal)])
print(round(obj_zal[!is.na(obj_zal)], 4))


