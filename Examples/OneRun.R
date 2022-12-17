set.seed(1)
rm(list = ls())
library(RColorBrewer)
source("TruncatedBeta/functions.R")


########################################
###### Example 1: Table-2 ##############
########################################

load(file = "Quadratic/Out/quad-objects_sq1e3.Rdata")

print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

print(round(quantile(time_bqn1, c(.5, .25, .75)), 3))
print(quantile(eval_bqn1, c(.5, .25, .75)))
print(round(quantile(obj_bqn1, c(.5, .25, .75)), 4))

print(round(quantile(time_bqn2, c(.5, .25, .75)), 3))
print(quantile(eval_bqn2, c(.5, .25, .75)))
print(round(quantile(obj_bqn2, c(.5, .25, .75)), 4))

print(round(quantile(time_bqn3, c(.5, .25, .75)), 3))
print(quantile(eval_bqn3, c(.5, .25, .75)))
print(round(quantile(obj_bqn3, c(.5, .25, .75)), 4))

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

print(round(quantile(time_zal2, c(.5, .25, .75)), 3))
print(quantile(eval_zal2, c(.5, .25, .75)))
print(round(quantile(obj_zal2, c(.5, .25, .75)), 4))

print(round(quantile(time_dar, c(.5, .25, .75)), 3))
print(quantile(eval_dar, c(.5, .25, .75)))
print(round(quantile(obj_dar, c(.5, .25, .75)), 4))

########################################
###### Example 2: Table-3 #############
########################################

load(file = "TruncatedBeta/Out/beta-objects1.Rdata")

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

########################################
###### Example 2: Figure-3 #############
########################################

load(file = "TruncatedBeta/Out/beta-objects1.Rdata")

batch <- 4
freq1 <- 15
freq2 <- 5
freq3 <- 2
freq4 <- 2
x <- seq(0.000001, 0.6, .01)
y <- seq(0, 2, 0.01)
z <- outer(X=x, Y=y, f, batch=4, freq1=freq1, freq2=freq2, freq3=freq3, freq4=freq4)


pdf(file = "TruncatedBeta/Out/beta-contour1_MM.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain_mm[1:iter_mm,1],chain_mm[1:iter_mm,2], col = c(rep(1,(iter_mm-1)), 2), pch = c(rep(1,(iter_mm-1)), 19), cex = c(rep(2,(iter_mm-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_BQN1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn1$p.inter[,1],fp_bqn1$p.inter[,2], col = c(rep(1,(fp_bqn1$iter-1)), 2), pch = c(rep(1,(fp_bqn1$iter-1)), 19), cex = c(rep(2,(fp_bqn1$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_BQN2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_bqn2$p.inter[,1],fp_bqn2$p.inter[,2], col = c(rep(1,(fp_bqn2$iter-1)), 2), pch = c(rep(1,(fp_bqn2$iter-1)), 19), cex = c(rep(2,(fp_bqn2$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_LBQN.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_lbqn$p.inter[,1],fp_lbqn$p.inter[,2], col = c(rep(1,(fp_lbqn$iter-1)), 2), pch = c(rep(1,(fp_lbqn$iter-1)), 19), cex = c(rep(2,(fp_lbqn$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_SqS1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2);
  points(fp_sq1_img$p.inter[,1],fp_sq1_img$p.inter[,2], col = c(rep(1,(fp_sq1_img$iter-1)), 2), pch = c(rep(1,(fp_sq1_img$iter-1)), 19), cex = c(rep(2,(fp_sq1_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_SqS2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq2_img$p.inter[,1],fp_sq2_img$p.inter[,2], col = c(rep(1,(fp_sq2_img$iter-1)), 2), pch = c(rep(1,(fp_sq2_img$iter-1)), 19), cex = c(rep(2,(fp_sq2_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_SqS3.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_sq3_img$p.inter[,1],fp_sq3_img$p.inter[,2], col = c(rep(1,(fp_sq3_img$iter-1)), 2), pch = c(rep(1,(fp_sq3_img$iter-1)), 19), cex = c(rep(2,(fp_sq3_img$iter-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_ZAL.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal_temp$Xhist[1,],fp_zal_temp$Xhist[2,], col = c(rep(1,(ncol(fp_zal_temp$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal_temp$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal_temp$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_ZAL2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_zal_temp2$Xhist[1,],fp_zal_temp2$Xhist[2,], col = c(rep(1,(ncol(fp_zal_temp2$Xhist)-1)), 2), pch = c(rep(1,(ncol(fp_zal_temp2$Xhist)-1)), 19), cex = c(rep(2,(ncol(fp_zal_temp2$Xhist)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

pdf(file = "TruncatedBeta/Out/beta-contour1_DAR.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp_dar$p.intermed[,2],fp_dar$p.intermed[,3], col = c(rep(1,(nrow(fp_dar$p.intermed)-1)), 2), pch = c(rep(1,(ncol(fp_dar$p.intermed)-1)), 19), cex = c(rep(2,(ncol(fp_dar$p.intermed)-1)), 2.5))}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE), xlab = expression(pi), ylab = expression(alpha))
dev.off()

########################################
###### Example 3: Table-4 #############3ewsaz
########################################

D = c("descent", "ascent")

for (d in 1:2)
{
  dir <- D[d]
  load(paste("GenEigen/Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))

  print(paste("Case = ", dir))
  print(paste("MM Algo ---------- Time: ", round(quantile(time_mm[1], .5), 3),  "Iterations: ", quantile(eval_mm[1], .5), "Negative log likelihood: ", round(quantile(obj_mm[1], .5), 4)))
  print(paste("BQN1 ------------- Time: ", round(quantile(time_bqn1[1], .5), 3), "Iterations: ", quantile(eval_bqn1[1], .5), "Negative log likelihood: ", round(quantile(obj_bqn1[1], .5), 4)))
  print(paste("BQN2 ------------- Time: ", round(quantile(time_bqn2[1], .5), 3), "Iterations: ", quantile(eval_bqn2[1], .5), "Negative log likelihood: ", round(quantile(obj_bqn2[1], .5), 4)))
  print(paste("BQN3 ------------- Time: ", round(quantile(time_bqn3[1], .5), 3), "Iterations: ", quantile(eval_bqn3[1], .5), "Negative log likelihood: ", round(quantile(obj_bqn3[1], .5), 4)))
  print(paste("LBQN ------------- Time: ", round(quantile(time_lbqn[1], .5), 3), "Iterations: ", quantile(eval_lbqn[1], .5), "Negative log likelihood: ", round(quantile(obj_lbqn[1], .5), 4)))
  print(paste("SQUAREM1 --------- Time: ", round(quantile(time_sq1[1], .5), 3), "Iterations: ", quantile(eval_sq1[1], .5), "Negative log likelihood: ", round(quantile(obj_sq1[1], .5), 4)))
  print(paste("SQUAREM2 --------- Time: ", round(quantile(time_sq2[1], .5), 3), "Iterations: ", quantile(eval_sq2[1], .5), "Negative log likelihood: ", round(quantile(obj_sq2[1], .5), 4)))
  print(paste("SQUAREM3 --------- Time: ", round(quantile(time_sq3[1], .5), 3), "Iterations: ", quantile(eval_sq3[1], .5), "Negative log likelihood: ", round(quantile(obj_sq3[1], .5), 4)))
  print(paste("ZAL, q=1 --------- Time: ", round(quantile(time_zal[1], .5), 3), "Iterations: ", quantile(eval_zal[1], .5), "Negative log likelihood: ", round(quantile(obj_zal[1], .5), 4)))
  print(paste("ZAL, q=2 --------- Time: ", round(quantile(time_zal2[1], .5), 3), "Iterations: ", quantile(eval_zal2[1], .5), "Negative log likelihood: ", round(quantile(obj_zal2[1], .5), 4)))
  print(paste("ZAL, q=min(p,10)-- Time: ", round(quantile(time_zal3[1], .5), 3), "Iterations: ", quantile(eval_zal3[1], .5), "Negative log likelihood: ", round(quantile(obj_zal3[1], .5), 4)))
  print(paste("DAAREM ----------- Time: ", round(quantile(time_dar[1], .5), 3), "Iterations: ", quantile(eval_dar[1], .5), "Negative log likelihood: ", round(quantile(obj_dar[1], .5), 4)))
}

########################################
###### Example 4: Table - 5 #############
########################################

load(file = "MultiT/Out/multiT-objects.Rdata")

print(paste("Number of failures:", sum(is.na(time_mm))))
print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_pxem))))
print(round(quantile(time_pxem[!is.na(time_pxem)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_pxem[!is.na(eval_pxem)], probs = c(.5, .25, .75)))
print(round(quantile(obj_pxem[!is.na(obj_pxem)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_bqn1))))
print(round(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_bqn2))))
print(round(quantile(time_bqn2[!is.na(time_bqn2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn2[!is.na(eval_bqn2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn2[!is.na(obj_bqn2)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_bqn3))))
print(round(quantile(time_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn3[!is.na(eval_bqn3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn3[!is.na(obj_bqn3)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_lbqn))))
print(round(quantile(time_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_lbqn[!is.na(eval_lbqn)], probs = c(.5, .25, .75)))
print(round(quantile(obj_lbqn[!is.na(obj_lbqn)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_sq1))))
print(round(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq1[!is.na(eval_sq1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq1[!is.na(obj_sq1)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_sq2))))
print(round(quantile(time_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq2[!is.na(eval_sq2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq2[!is.na(obj_sq2)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_sq3))))
print(round(quantile(time_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq3[!is.na(eval_sq3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq3[!is.na(obj_sq3)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_zal))))
print(round(quantile(time_zal[!is.na(time_zal)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal[!is.na(eval_zal)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal[!is.na(obj_zal)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_zal2))))
print(round(quantile(time_zal2[!is.na(time_zal2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal2[!is.na(eval_zal2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal2[!is.na(obj_zal2)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_zal3))))
print(round(quantile(time_zal3[!is.na(time_zal3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_zal3[!is.na(eval_zal3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_zal3[!is.na(obj_zal3)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(time_dar))))
print(round(quantile(time_dar, probs = c(.5, .25, .75)), 3))
print(quantile(eval_dar, probs = c(.5, .25, .75)))
print(round(quantile(obj_dar, probs = c(.5, .25, .75)), 4))

########################################
###### Example 2: Appendix:Table-2 #####
########################################

for (d in 2:4)
{
  load(file = paste("TruncatedBeta/Out/beta-objects", d, ".Rdata", sep=""))
  print(paste("Dataset:", d))
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

########################################
###### Example 5: Appendix:Table-4 #####
########################################


load(file = "Poisson/Out/poisson-objects.Rdata")

print(paste("Number of failures:", sum(is.na(eval_mm))))
print(round(quantile(time_mm, probs = c(.5, .25, .75)), 3))
print(quantile(eval_mm, probs = c(.5, .25, .75)))
print(round(quantile(obj_mm, probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_bqn1))))
print(round(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_bqn2))))
print(round(quantile(time_bqn2[!is.na(time_bqn2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_bqn2[!is.na(eval_bqn2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_bqn2[!is.na(obj_bqn2)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_lbqn))))
print(round(quantile(time_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_lbqn[!is.na(eval_lbqn)], probs = c(.5, .25, .75)))
print(round(quantile(obj_lbqn[!is.na(obj_lbqn)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_sq1))))
print(round(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq1[!is.na(eval_sq1)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq1[!is.na(obj_sq1)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_sq2))))
print(round(quantile(time_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq2[!is.na(eval_sq2)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq2[!is.na(obj_sq2)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_sq3))))
print(round(quantile(time_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)), 3))
print(quantile(eval_sq3[!is.na(eval_sq3)], probs = c(.5, .25, .75)))
print(round(quantile(obj_sq3[!is.na(obj_sq3)], probs = c(.5, .25, .75)), 4))

print(paste("Number of failures:", sum(is.na(eval_zal))))

print(paste("Number of failures:", sum(is.na(eval_zal2))))

print(paste("Number of failures:", sum(is.na(eval_dar))))
print(round(quantile(time_dar, probs = c(.5, .25, .75)), 3))
print(quantile(eval_dar, probs = c(.5, .25, .75)))
print(round(quantile(obj_dar, probs = c(.5, .25, .75)), 4))

