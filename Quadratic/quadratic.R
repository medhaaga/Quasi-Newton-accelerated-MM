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

dim <- 2
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

N <- 1
time_rep.mm <- rep(0,N)
fevals.mm <- rep(0,N)
obj.values.mm <- rep(0, N)
chain <- as.matrix(start, nrow = 1, ncol = dim)

for (j in 1:N)
{
  print(j)

  start <- truth + as.matrix(start.all[j,])
  now <- start
  new <- start
  iter <- 1
  diff <- 100


  start.time <- Sys.time()
  while(diff > tol)
  {
    if (iter %% 100 == 0) print(diff)
    new <- update(now, A, b, L)
    chain <- cbind(chain, new)
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
  }
  end.time <- Sys.time()
  pdf(file = "Out/quad-contour_MM.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain[1,],chain[2,])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()

  time_rep.mm [j] <- end.time - start.time
  fevals.mm[j] <- iter
  obj.values.mm[j] <- objective(new, A, b, L)
}

print(quantile(time_rep.mm, probs = c(0, .5, 1)))
print(quantile(fevals.mm, probs = c(0, .5, 1)))
print(quantile(obj.values.mm, probs = c(0, .5, 1)))

##################################################
#### BFGS
##################################################


N <- 1
time_rep.bfgs <- rep(0,N)
fevals.bfgs <- rep(0,N)
levals <- rep(0, N)
obj.values.bfgs <- rep(0, N)
fails <- 0

for (j in 1:N)
{
  print(j)
  start <- truth + as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- BFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(tol = tol, objfn.inc = .001, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_BFGS.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()

  if(fp$convergence){
    time_rep.bfgs [j] <- (end.time - start.time)*60
    fevals.bfgs[j] <- fp$fpevals
    obj.values.bfgs[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.bfgs [j] <- NA
    fevals.bfgs[j] <- NA
    obj.values.bfgs[j] <- NA
  }
}

print(quantile(time_rep.bfgs*60, probs = c(0, .5, 1)))
print(quantile(fevals.bfgs, probs = c(0, .5, 1)))
print(quantile(obj.values.bfgs, probs = c(0, .5, 1)))

#####################################################
#### L-BFGS
#####################################################


N <- 1
time_rep.lbfgs <- rep(0,N)
fevals.lbfgs <- rep(0,N)
levals <- rep(0, N)
obj.values.lbfgs <- rep(0, N)
tol = 1e-3
fails <- 0

for (j in 1:N)
{
  print(j)
  start <- truth + as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- LBFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(tol = tol, objfn.inc = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_LBFGS.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()
  if(fp$convergence){
    time_rep.lbfgs [j] <- (end.time - start.time)*60
    fevals.lbfgs[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.lbfgs[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.lbfgs [j] <- NA
    fevals.lbfgs[j] <- NA
    levals[j] <- NA
    obj.values.lbfgs[j] <- NA
  }
}

print(quantile(time_rep.lbfgs, probs = c(0, .5, 1)))
print(quantile(fevals.lbfgs, probs = c(0, .5, 1)))
print(quantile(obj.values.lbfgs, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 1
##########################################

N <- 1
time_rep.sq1 <- rep(0,N)
fevals.sq1 <- rep(0,N)
levals <- rep(0, N)


obj.values.sq1 <- rep(0, N)
fails <- 0


for (j in 1:N){
  print(j)
  start <- truth + as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_SqS1.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()
  if(fp$convergence){
    time_rep.sq1 [j] <- end.time - start.time
    fevals.sq1[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.sq1[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.sq1 [j] <- NA
    fevals.sq1[j] <- NA
    levals[j] <- NA
    obj.values.sq1[j] <- NA
  }
}
print(quantile(time_rep.sq1, probs = c(0, .5, 1)))
print(quantile(fevals.sq1, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.sq1, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 2
##########################################

N <- 1
time_rep.sq2 <- rep(0,N)
fevals.sq2 <- rep(0,N)
levals <- rep(0, N)
obj.values.sq2 <- rep(0, N)
tol = 1e-3
fails <- 0


for (j in 1:N){
  print(j)
  start <- truth + as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed= TRUE))
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_SqS2.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()

  if(fp$convergence){
    time_rep.sq2 [j] <- end.time - start.time
    fevals.sq2[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.sq2[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.sq2 [j] <- NA
    fevals.sq2[j] <- NA
    levals[j] <- NA
    obj.values.sq2[j] <- NA
  }
}
print(quantile(time_rep.sq2, probs = c(0, .5, 1)))
print(quantile(fevals.sq2, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.sq2, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 3
##########################################

N <- 1
time_rep.sq3 <- rep(0,N)
fevals.sq3 <- rep(0,N)
levals.sq3 <- rep(0, N)
obj.values.sq3 <- rep(0, N)
fails <- 0


for (j in 1:N){
  print(j)
  start <-truth +as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- squarem(start, A=A, b=b, L=L, fixptfn = update, objfn = objective, control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed= TRUE))
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_SqS3.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()
  if(fp$convergence){
    time_rep.sq3 [j] <- end.time - start.time
    fevals.sq3[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.sq3[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.sq3 [j] <- NA
    fevals.sq3[j] <- NA
    levals[j] <- NA
    obj.values.sq3[j] <- NA
  }
}
print(quantile(time_rep.sq3, probs = c(0, .5, 1)))
print(quantile(fevals.sq3, probs = c(0, .5, 1)))
print(quantile(levals.sq3, probs = c(0, .5, 1)))
print(quantile(obj.values.sq3, probs = c(0, .5, 1)))

#################################################
##### ZAL; q=2
#################################################

q <- 2
N <- 1
time_rep.zal <- rep(0,N)
fevals.zal <- rep(0,N)
levals <- rep(0,N)
obj.values.zal <- rep(0, N)
accept <- rep(0, N)
reject <- rep(0, N)
fails <- 0


for (j in 1:N){
  print(j)
  start <- truth + as.matrix(start.all[j,])

  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, qn=q, fx_obj = objective, max_iter = 5e4, tol = tol, A=A, b=b, L=L)
  end.time <- Sys.time()

  pdf(file = "Out/quad-contour_ZAL.pdf", height = 5, width = 7)
  filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$Xhist[1,],fp$Xhist[2,])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
  dev.off()
  if (fp$convergence == TRUE){
    time_rep.zal[j] <- end.time - start.time
    fevals.zal[j] <- fp$fevals
    levals[j] <- fp$levals
    obj.values.zal[j] <- fp$objective
    accept[j] <- fp$accept
    reject[j] <- fp$reject}
  else{
    time_rep.zal[j] <- NA
    fevals.zal[j] <- NA
    levals[j] <- NA
    obj.values.zal[j] <- NA
    accept[j] <- NA
    reject[j] <- NA
  }
}

print(quantile(time_rep.zal, probs = c(0, .5, 1)))
print(quantile(fevals.zal, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.zal, probs = c(0, .5, 1)))
print(quantile(accept, probs = c(0, .5, 1)))
print(quantile(reject, probs = c(0, .5, 1)))


####################################
#### Boxplots and scatterplots
#####################################

save(obj.values.bfgs, obj.values.mm, obj.values.sq1, obj.values.sq2,
     obj.values.sq3, obj.values.zal, fevals.mm, fevals.bfgs, fevals.sq1,
     fevals.sq2, fevals.sq3, fevals.zal, time_rep.mm, time_rep.bfgs,
     time_rep.sq1, time_rep.sq2, time_rep.sq3, time_rep.zal, file = "Out/quad_sd100.Rdata")

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
