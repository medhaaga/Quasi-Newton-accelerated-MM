set.seed(1)
library(pracma)
library(SQUAREM)
library(BfgsQN)
source("qnamm.r")
library(RColorBrewer)

objective <- function(x, A, b, L)
{
  obj <- 0.5*crossprod(x, (A %*% x)) + crossprod(b, x)
  return (obj)
}

update <- function(now, A, b, L)
{
  new <- now - ((A %*% now) + b)/L
  return (new)
}

f <- function(x, y, A, b)
{
  return ((A[1,1]*x^2) + x*y*(A[1,2]+A[2,1]) + A[2,2]*y^2 + b[1]*x + b[2]*y)
}

#######################################

set.seed(1)
dim <- 5
u <- matrix(rnorm(dim*dim, mean = 0, sd = 10), dim, dim)
A <- t(u) %*% u
L <- 1.001*max(abs(eigen(A)$values))
b <- as.matrix(rnorm(dim, mean = 0, sd = 10))
truth <- - solve(A) %*% b
true.obj <- objective(truth, A, b, L)
N <- 1
start.all <- matrix(rnorm(dim*N, mean = 0, sd = 10), nrow = N, ncol = dim)
tol = 1e-7

x.range <- abs(truth[1] - start.all[1,1])
y.range <- abs(truth[2] - start.all[1,2])
x <- seq(-5, truth[2] + x.range + 3, .1)
y <- seq(truth[2] - y.range - 3, 1, .1)
z <- outer(X=x, Y=y, f, A=A, b=b)


############################################

###############################################
#### MM Algorithm
###############################################

time.mm <- rep(0,N)
eval.mm <- rep(0,N)
obj.mm <- rep(0, N)

for (j in 1:N)
{
  print(j)

  start <- truth + as.matrix(start.all[j,])
  chain <- as.matrix(t(start), nrow = 1, ncol = dim)
  
  now <- start
  new <- start
  iter <- 1
  diff <- 100


  start.time <- Sys.time()
  while(diff > tol)
  {
    if (iter %% 100 == 0) print(diff)
    new <- update(now, A, b, L)
    chain <- rbind(chain, t(new))
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
  }
  end.time <- Sys.time()

  time.mm [j] <- end.time - start.time
  eval.mm[j] <- iter
  obj.mm[j] <- objective(new, A, b, L)
}

#pdf(file = "Out/quad-contour_MM.pdf", height = 5, width = 7)
#filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain[,1],chain[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
#dev.off()

print(quantile(time.mm, probs = c(0, .5, 1)))
print(quantile(eval.mm, probs = c(0, .5, 1)))
print(quantile(obj.mm, probs = c(0, .5, 1)))

##################################################
#### BFGS
##################################################

time_bfgs <- rep(0, N)
obj_bfgs <- rep(0, N)
eval_bfgs <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- BFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(tol = tol, obj.tol=1e-7, objfn.inc = .001, step.max=100, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_bfgs[i] <- end.time - start.time
  obj_bfgs[i] <- fp$value.objfn
  eval_bfgs[i] <- fp$fpevals

}

pdf(file = "Out/quad-contour_BFGS.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(quantile(time_bfgs, probs = c(0, .5, 1)))
print(quantile(eval_bfgs, probs = c(0, .5, 1)))
print(quantile(obj_bfgs, probs = c(0, .5, 1)))

#####################################################
#### L-BFGS
#####################################################


time_lbfgs <- rep(0, N)
obj_lbfgs <- rep(0, N)
eval_lbfgs <- rep(0, N)


for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- LBFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(tol = tol, obj.tol=1e-7, objfn.inc = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()


  time_lbfgs[i] <- end.time - start.time
  obj_lbfgs[i] <- fp$value.objfn
  eval_lbfgs[i] <- fp$fpevals

}
pdf(file = "Out/quad-contour_LBFGS.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(quantile(time_lbfgs, probs = c(0, .5, 1)))
print(quantile(eval_lbfgs, probs = c(0, .5, 1)))
print(quantile(obj_lbfgs, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 1
##########################################

eval_sq1 <- rep(0, N)
time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()


  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
}

pdf(file = "Out/quad-contour_SqS1.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()

print(quantile(time_sq1, probs = c(0, .5, 1)))
print(quantile(eval_sq1, probs = c(0, .5, 1)))
print(quantile(obj_sq1, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 2
##########################################

eval_sq2 <- rep(0, N)
time_sq2 <- rep(0, N)
obj_sq2 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()


  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
}

pdf(file = "Out/quad-contour_SqS2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(quantile(time_sq2, probs = c(0, .5, 1)))
print(quantile(eval_sq2, probs = c(0, .5, 1)))
print(quantile(obj_sq2, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 3
##########################################

eval_sq3 <- rep(0, N)
time_sq3 <- rep(0, N)
obj_sq3 <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective,
                control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()


  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
}

pdf(file = "Out/quad-contour_SqS2.pdf", height = 5, width = 7)
filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
dev.off()
print(quantile(time_sq3, probs = c(0, .5, 1)))
print(quantile(eval_sq3, probs = c(0, .5, 1)))
print(quantile(obj_sq3, probs = c(0, .5, 1)))

#################################################
##### ZAL; q=2
#################################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- start.all[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = objective, qn=2, A=A, b=b, L=L, max_iter = 1e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
}

fails <- sum(is.na(obj_zal))
time_zal <- time_zal[!is.na(time_zal)]
eval_zal <- eval_zal[!is.na(eval_zal)]
obj_zal <- obj_zal[!is.na(obj_zal)]
print(quantile(time_zal, probs = c(0, .5, 1)))
print(quantile(eval_zal, probs = c(0, .5, 1)))
print(quantile(obj_zal, probs = c(0, .5, 1)))


####################################
#### Boxplots and scatterplots
#####################################

save(time_mm, time_bfgs, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bfgs, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bfgs, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = "Out/objects.Rdata")

time_range <- range(time_mm, time_bfgs, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal)
eval_range <- range(eval_mm, eval_bfgs, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal)
obj_range <- range(obj_mm, obj_bfgs, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal)

df1 <- data.frame("MM" = obj.values.mm, "SqS1" = obj.values.sq1, "SqS2" = obj.values.sq2,
                 "SqS3" = obj.values.sq3, "BFGS" = obj.values.bfgs, "ZAL" = obj.values.zal)

df2 <- data.frame("MM" = fevals.mm, "SqS1" = fevals.sq1, "SqS2" = fevals.sq2,
                  "SqS3" = fevals.sq3, "BFGS" = fevals.bfgs, "ZAL" = fevals.zal)

df3 <- data.frame("MM" = time_rep.mm, "SqS1" = time_rep.sq1, "SqS2" = time_rep.sq2,
                  "SqS3" = time_rep.sq3, "BFGS" = time_rep.bfgs, "ZAL" = time_rep.zal)


load(file = "Out/quad_sd1.Rdata")

pdf(file = "Out/boxplot_quad_s100.pdf", width = 12, height = 6)
par(mfrow = c(1,3))
boxplot(df1, xlab = "Algorithm", ylab = "Objective value")
boxplot(df2, xlab = "Algorithm", ylab = "No. of evaluations")
boxplot(df3, xlab = "Algorithm", ylab = "Time")
dev.off()

pdf(file = "Out/scatter_plot_quad_s1.pdf")
plot(fevals.mm, obj.values.mm, pch = 19, col = "cadetblue1", xlab = "Number of updates", ylab = "Objective Value",
     xlim = range(fevals.mm, fevals.bfgs, fevals.lbfgs,fevals.sq1, fevals.sq2, fevals.sq3, fevals.zal),
     ylim = range(obj.values.bfgs, obj.values.lbfgs, obj.values.mm, obj.values.sq1, obj.values.sq2, obj.values.sq3, obj.values.zal))
points(fevals.bfgs, obj.values.bfgs, pch = 19, col = "cyan3")
points(fevals.lbfgs, obj.values.lbfgs, pch = 19, col = "hotpink")
points(fevals.sq1, obj.values.sq1, pch = 19, col = "lightpink")
points(fevals.sq2, obj.values.sq2, pch = 19, col = "coral")
points(fevals.sq3, obj.values.sq3, pch = 19, col = "blue")
points(fevals.zal, obj.values.zal, pch = 19, col = "mediumorchid4")
legend("bottomright", legend = c("MM", "BFGS", "L-BFGS", "SqS1", "SqS2", "SqS3", "ZAL"), col = c("cadetblue1", "cyan3", "hotpink", "lightpink", "coral", "blue", "mediumorchid4"), pch = 19)
dev.off()
