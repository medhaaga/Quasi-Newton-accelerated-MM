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

set.seed(10)
dim <- 5
u <- matrix(rnorm(dim*dim, mean = 0, sd = 10), dim, dim)
A <- t(u) %*% u
L <- 1.001*max(abs(eigen(A)$values))
b <- as.matrix(rnorm(dim, mean = 0, sd = 10))
truth <- - solve(A) %*% b
true.obj <- objective(truth, A, b, L)
N <- 10
start.all <- matrix(rnorm(dim*N, mean = 0, sd = 1e3), nrow = N, ncol = dim)
tol = 1e-7

#x.range <- abs(truth[1] - start.all[1,1])
#y.range <- abs(truth[2] - start.all[1,2])
#x <- seq(-5, truth[2] + x.range + 3, .1)
#y <- seq(truth[2] - y.range - 3, 1, .1)
#z <- outer(X=x, Y=y, f, A=A, b=b)


############################################

###############################################
#### MM Algorithm
###############################################

time_mm <- rep(0,N)
eval_mm <- rep(0,N)
obj_mm <- rep(0, N)

for (j in 1:N)
{
  print(j)

  start <- truth + as.matrix(start.all[j,])
  objval <- objective(start, A=A, b=b, L=L)
  chain <- as.matrix(t(start), nrow = 1, ncol = dim)
  chain_mm <- objval
  now <- start
  new <- start
  iter <- 1
  diff <- 100


  start.time <- Sys.time()
  while(diff > tol)
  {
    if (iter %% 1000 == 0) print(diff)
    new <- update(now, A, b, L)
    objval <- objective(new, A=A, b=b, L=L)
    chain_mm <- c(chain_mm, objval)
    chain <- rbind(chain, t(new))
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
  }
  end.time <- Sys.time()

  time_mm [j] <- end.time - start.time
  eval_mm[j] <- iter
  obj_mm[j] <- objective(new, A, b, L)
  
}

#pdf(file = "Out/quad-contour_MM.pdf", height = 5, width = 7)
#filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(chain[,1],chain[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
#dev.off()

print(quantile(time_mm, probs = c(.1, .5, .9)))
print(quantile(eval_mm, probs = c(.1, .5, .9)))
print(quantile(obj_mm, probs = c(.1, .5, .9)))

##################################################
#### BFGS, q=1
##################################################

time_bfgs1 <- rep(0, N)
obj_bfgs1 <- rep(0, N)
eval_bfgs1 <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- BFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, 
             control = list(qn=1, tol=tol, objfn.inc=100, step.max=100000, step.min=1, maxiter=5e4, intermed=TRUE))
  end.time <- Sys.time()

  time_bfgs1[i] <- end.time - start.time
  obj_bfgs1[i] <- fp$value.objfn
  eval_bfgs1[i] <- fp$fpevals
  chain_bfgs1 <- fp$p.inter[,(dim+1)]

}

#pdf(file = "Out/quad-running_BFGS.pdf", height = 5, width = 7)
#plot(20:fp$iter, fp$p.inter[20:fp$iter,11], type = "l")
#dev.off()

print(quantile(time_bfgs1, probs = c(.1, .5, .9)))
print(quantile(eval_bfgs1, probs = c(.1, .5, .9)))
print(quantile(obj_bfgs1, probs = c(.1, .5, .9)))

##################################################
#### BFGS, q=2
##################################################

time_bfgs2 <- rep(0, N)
obj_bfgs2 <- rep(0, N)
eval_bfgs2 <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])
  
  start.time <- Sys.time()
  fp <- BFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, 
             control = list(qn=2, tol=tol, objfn.inc=1, step.max=1000, step.min=1, maxiter=5e4, intermed=TRUE))
  end.time <- Sys.time()
  
  time_bfgs2[i] <- end.time - start.time
  obj_bfgs2[i] <- fp$value.objfn
  eval_bfgs2[i] <- fp$fpevals
  chain_bfgs2 <- fp$p.inter[,(dim+1)]
  
}

#pdf(file = "Out/quad-running_BFGS.pdf", height = 5, width = 7)
#plot(20:fp$iter, fp$p.inter[20:fp$iter,11], type = "l")
#dev.off()

print(quantile(time_bfgs2, probs = c(.1, .5, .9)))
print(quantile(eval_bfgs2, probs = c(.1, .5, .9)))
print(quantile(obj_bfgs2, probs = c(.1, .5, .9)))

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
  fp <- LBFGS(par = start, A=A, b=b, L=L, fixptfn = update, objfn = objective, 
              control = list(tol = tol, objfn.inc = .01, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()


  time_lbfgs[i] <- end.time - start.time
  obj_lbfgs[i] <- fp$value.objfn
  eval_lbfgs[i] <- fp$fpevals
  chain_lbfgs <- fp$p.inter[,(dim+1)]
}
#pdf(file = "Out/quad-running_LBFGS.pdf", height = 5, width = 7)
#plot(20:fp$iter, fp$p.inter[20:fp$iter,11], type = "l")
#dev.off()
print(quantile(time_lbfgs, probs = c(.1, .5, .9)))
print(quantile(eval_lbfgs, probs = c(.1, .5, .9)))
print(quantile(obj_lbfgs, probs = c(.1, .5, .9)))

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
  chain_sq1 <- fp$p.inter[,(dim+1)]
}

#pdf(file = "Out/quad-contour_SqS1.pdf", height = 5, width = 7)
#filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
#dev.off()

print(quantile(time_sq1, probs = c(.1, .5, .9)))
print(quantile(eval_sq1, probs = c(.1, .5, .9)))
print(quantile(obj_sq1, probs = c(.1, .5, .9)))

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
  chain_sq2 <- fp$p.inter[,(dim+1)]
}

#pdf(file = "Out/quad-contour_SqS2.pdf", height = 5, width = 7)
#filled.contour(x,y,z,plot.axes = { axis(1); axis(2); points(fp$p.inter[,1],fp$p.inter[,2])}, color.palette = function(n) hcl.colors(n, "RdPu", rev = TRUE))
#dev.off()
print(quantile(time_sq2, probs = c(.1, .5, .9)))
print(quantile(eval_sq2, probs = c(.1, .5, .9)))
print(quantile(obj_sq2, probs = c(.1, .5, .9)))

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
  chain_sq3 <- fp$p.inter[,(dim+1)]
}

print(quantile(time_sq3, probs = c(.1, .5, .9)))
print(quantile(eval_sq3, probs = c(.1, .5, .9)))
print(quantile(obj_sq3, probs = c(.1, .5, .9)))

#################################################
##### ZAL; q=2
#################################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- start.all[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = objective, qn=2, A=A, b=b, L=L, max_iter = 1e5, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
  chain_zal <- fp$Xhist[dim+1,]
}

fails <- sum(is.na(obj_zal))
time_zal <- time_zal[!is.na(time_zal)]
eval_zal <- eval_zal[!is.na(eval_zal)]
obj_zal <- obj_zal[!is.na(obj_zal)]
print(quantile(time_zal, probs = c(.1, .5, .9)))
print(quantile(eval_zal, probs = c(.1, .5, .9)))
print(quantile(obj_zal, probs = c(.1, .5, .9)))


####################################
#### Boxplots and scatterplots
#####################################

save(time_mm, time_bfgs1, time_bfgs2, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bfgs1, eval_bfgs2, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bfgs1, obj_bfgs2, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal, file = "Out/quad-objects_sq1e3.Rdata")

time_range <- range(time_mm, time_bfgs1, time_bfgs2, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal)
eval_range <- range(eval_mm, eval_bfgs1, eval_bfgs2, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal)
obj_range <- range(obj_mm, obj_bfgs1, obj_bfgs2, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal)

df1 <- data.frame("BFGS1" = obj_bfgs1, "BFGS2" = obj_bfgs2, "LBFGS" = obj_lbfgs, 
                  "Sq1" = obj_sq1, "Sq2" = obj_sq2, "Sq3" = obj_sq3, "ZAL" = obj_zal)

df2 <- data.frame("BFGS1" = eval_bfgs1, "BFGS2" = eval_bfgs2, "LBFGS" = eval_lbfgs, 
                  "Sq1" = eval_sq1, "Sq2" = eval_sq2, "Sq3" = eval_sq3, "ZAL" = eval_zal)

df3 <- data.frame("BFGS1" = time_bfgs1, "BFGS2" = time_bfgs2, "LBFGS" = time_lbfgs, 
                  "Sq1" = time_sq1, "Sq2" = time_sq2, "Sq3" = time_sq3, "ZAL" = time_zal)


pdf(file = "Out/quad-boxplot_sd1e3.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
boxplot(df2, xlab = "Algorithm", ylab = "No. of evaluations")
boxplot(df3, xlab = "Algorithm", ylab = "Time")
dev.off()
