##################################################
########## Quadratic Minimization ################
##################################################
set.seed(1)
rm(list = ls())
library(pracma)
library(quasiNewtonMM)
library(turboEM)
library(daarem)
library(RColorBrewer)

objective <- function(x, A, a, L)
{
  obj <- 0.5*crossprod(x, (A %*% x)) + crossprod(a, x)
  return (obj)
}

neg.objective <- function(x, A, a, L){
  obj <- 0.5*crossprod(x, (A %*% x)) + crossprod(a, x)
  return (-obj)
}

update <- function(now, A, a, L)
{
  new <- now - ((A %*% now) + a)/L
  return (new)
}

#######################################


dim <- 100
u <- matrix(rnorm(dim*dim, mean = 0, sd = 10), dim, dim)
A <- t(u) %*% u
L <- 1.001*norm(max(abs(eigen(A)$values)), type = "2")
a <- as.matrix(rnorm(dim, mean = 0, sd = 10))
truth <- - solve(A) %*% a
true.obj <- objective(truth, A, a, L)
N <- 10
start.all <- matrix(rnorm(dim*N, mean = 0, sd = 1000), nrow = N, ncol = dim)
tol = 1e-7


###############################################
#### Unaccelerated MM Algorithm
###############################################

time_mm <- rep(0,N)
eval_mm <- rep(0,N)
obj_mm <- rep(0, N)

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
    new <- update(now, A, a, L)
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
  }
  end.time <- Sys.time()

  time_mm [j] <- end.time - start.time
  eval_mm[j] <- iter
  obj_mm[j] <- objective(new, A, a, L)

}

print(round(quantile(time_mm, c(.5, .25, .75)), 3))
print(quantile(eval_mm, c(.5, .25, .75)))
print(round(quantile(obj_mm, c(.5, .25, .75)), 5))

#########################################
#### J & J QN1, Broyden's good method ####
#########################################

broy_fun <- function(x, func, ...){
  return(func(x, ...) - x)
}

time_qn1 <- rep(NA, N)
obj_qn1 <- rep(NA, N)
eval_qn1 <- rep(NA, N)

for (j in 1:N){
  print(j)
  start <- truth + as.matrix(start.all[j,])
  now <- start
  new <- start
  G_now <- broy_fun(now, func = update, A=A, a=a, L=L)
  G_new <- G_now
  H <- -diag(dim)
  itr <- 1
  diff <- 100
  start.time <- Sys.time()
  while(diff > tol){
    itr <- itr+1
    new <- now - H%*%G_now
    l_new <- objective(new, A, a, L)
    if(is.na(l_new)){
      print("Falling back to MM step")
      new <- G_now + now}
    G_new <- broy_fun(new, func = update, A=A, a=a, L=L)
    foo <- H%*%(G_new - G_now)
    H <- H + (((new-now) - foo)/as.numeric(t(new-now) %*% foo))%*%(t(new-now)%*%H)
    diff <- norm(new-now, "2")
    now <- new
    G_now <- G_new
  }
  end.time <- Sys.time()

  time_qn1[j] <- end.time - start.time
  obj_qn1[j] <- objective(new, A=A, a=a, L=L)
  eval_qn1[j] <- itr
}

print(quantile(time_qn1, probs = c(.5, .25, .75)))
print(quantile(eval_qn1, probs = c(.5, .25, .75)))
print(quantile(obj_qn1, probs = c(.5, .25, .75)))

##################################################
#### BQN, q=1
##################################################

time_bqn1 <- rep(0, N)
obj_bqn1 <- rep(0, N)
eval_bqn1 <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- BQN(par = start, A=A, a=a, L=L, fixptfn = update, objfn = objective,
             control = list(qn=1, tol=tol, maxiter=1e5))
  end.time <- Sys.time()

  time_bqn1[i] <- end.time - start.time
  obj_bqn1[i] <- fp$value.objfn
  eval_bqn1[i] <- fp$fpevals
}


print(round(quantile(time_bqn1, c(.5, .25, .75)), digits=3))
print(quantile(eval_bqn1, c(.5, .25, .75)))
print(round(quantile(obj_bqn1, c(.5, .25, .75)), 5))
# 
# ##################################################
#### BQN, q=2
##################################################

time_bqn2 <- rep(0, N)
obj_bqn2 <- rep(0, N)
eval_bqn2 <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- BQN(par = start, A=A, a=a, L=L, fixptfn = update, objfn = objective,
            control = list(qn=2, tol=tol, maxiter=1e5))
  end.time <- Sys.time()

  time_bqn2[i] <- end.time - start.time
  obj_bqn2[i] <- fp$value.objfn
  eval_bqn2[i] <- fp$fpevals
}


print(round(quantile(time_bqn2, c(.5, .25, .75)), digits=2))
print(quantile(eval_bqn2, c(.5, .25, .75)))
print(round(quantile(obj_bqn2, c(.5, .25, .75)), 5))

# ##################################################
# #### BQN, q=min(dim, 10)
# ##################################################

time_bqn3 <- rep(0, N)
obj_bqn3 <- rep(0, N)
eval_bqn3 <- rep(0, N)

for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- BQN(par = start, A=A, a=a, L=L, fixptfn = update, objfn = objective,
             control = list( qn=10, tol=1e-7, maxiter=1e5))
  end.time <- Sys.time()

  time_bqn3[i] <- end.time - start.time
  obj_bqn3[i] <- fp$value.objfn
  eval_bqn3[i] <- fp$fpevals
  }

print(round(quantile(time_bqn3, c(.5, .25, .75)), 3))
print(quantile(eval_bqn3, c(.5, .25, .75)))
print(round(quantile(obj_bqn3, c(.5, .25, .75)), 5))

#####################################################
#### L-BQN
#####################################################


time_lbqn <- rep(0, N)
obj_lbqn <- rep(0, N)
eval_lbqn <- rep(0, N)


for (i in 1:N)
{
  print(i)
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- LBQN(par = start, A=A, a=a, L=L, fixptfn = update, objfn = objective,
              control = list(tol = tol, maxiter = 1e5, m = 2, verbose=TRUE))
  end.time <- Sys.time()


  time_lbqn[i] <- end.time - start.time
  obj_lbqn[i] <- fp$value.objfn
  eval_lbqn[i] <- fp$fpevals
}

print(round(quantile(time_lbqn, c(.5, .25, .75)), 3))
print(quantile(eval_lbqn, c(.5, .25, .75)))
print(round(quantile(obj_lbqn, c(.5, .25, .75)), 5))

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
  fp <- fp <- turboem(par = start, fixptfn = update, objfn = objective, A=A, a=a, L=L, method = "squarem", control.method = list(K=1, version=1), control.run = list(tol=tol, maxiter=1e5))
  end.time <- Sys.time()


  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpeval
}

print(round(quantile(time_sq1, c(.5, .25, .75)), 3))
print(quantile(eval_sq1, c(.5, .25, .75)))
print(round(quantile(obj_sq1, c(.5, .25, .75)), 4))

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
  fp <- fp <- turboem(par = start, fixptfn = update, objfn = objective, A=A, a=a, L=L, method = "squarem", control.method = list(K=1, version=2), control.run = list(tol=tol, maxiter=1e5))
  end.time <- Sys.time()


  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpeval
}

print(round(quantile(time_sq2, c(.5, .25, .75)), 3))
print(quantile(eval_sq2, c(.5, .25, .75)))
print(round(quantile(obj_sq2, c(.5, .25, .75)), 4))

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
  fp <- turboem(par = start, fixptfn = update, objfn = objective, A=A, a=a, L=L, method = "squarem", control.method = list(K=1, version=3), control.run = list(tol=tol, maxiter=1e5))
  end.time <- Sys.time()


  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpeval
}

print(round(quantile(time_sq3, c(.5, .25, .75)), 3))
print(quantile(eval_sq3, c(.5, .25, .75)))
print(round(quantile(obj_sq3, c(.5, .25, .75)), 4))

#################################################
##### ZAL; q=2
#################################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = objective, A=A, a=a, L=L, method = "qn", control.method = list(qn=2), control.run = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  print(fp$convergence)
  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$value.objfn
  eval_zal[i] <- fp$fpeval
}

print(round(quantile(time_zal, c(.5, .25, .75)), 3))
print(quantile(eval_zal, c(.5, .25, .75)))
print(round(quantile(obj_zal, c(.5, .25, .75)), 4))

#################################################
##### ZAL; q=min(p,10)
#################################################

time_zal2 <- rep(NA, N)
obj_zal2 <- rep(NA, N)
eval_zal2 <- rep(NA, N)

for (i in 1:N){
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- turboem(par = start, fixptfn = update, objfn = objective, A=A, a=a, L=L, method = "qn", control.method = list(qn=min(dim,10)), control.run = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  print(fp$convergence)
  time_zal2[i] <- end.time - start.time
  obj_zal2[i] <- fp$value.objfn
  eval_zal2[i] <- fp$fpeval
}

print(round(quantile(time_zal2, c(.5, .25, .75)), 3))
print(quantile(eval_zal2, c(.5, .25, .75)))
print(round(quantile(obj_zal2, c(.5, .25, .75)), 4))

#################################################
##### DAAREM
#################################################

time_dar <- rep(NA, N)
obj_dar <- rep(NA, N)
eval_dar <- rep(NA, N)

for (i in 1:N){
  print(paste("Iter: ", i))
  start <- truth + as.matrix(start.all[i,])

  start.time <- Sys.time()
  fp <- daarem(par = start, fixptfn = update, objfn = neg.objective, A=A, a=a, L=L, control = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()
  time_dar[i] <- end.time - start.time
  obj_dar[i] <- -fp$value.objfn
  eval_dar[i] <- fp$fpeval
}

print(round(quantile(time_dar, c(.5, .25, .75)), 3))
print(quantile(eval_dar, c(.5, .25, .75)))
print(round(quantile(obj_dar, c(.5, .25, .75)), 4))


save(time_mm, time_bqn1, time_bqn2, time_bqn3, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, time_zal2, time_dar,
     eval_mm, eval_bqn1, eval_bqn2, eval_bqn3, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal, eval_zal2, eval_dar,
     obj_mm, obj_bqn1, obj_bqn2, obj_bqn3, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, obj_zal2, obj_dar, file = paste("Out/quad-objects_sq1e3.Rdata"))

####################################
#### Boxplot
#####################################


load(file = "Out/quad-objects_sq1e3.Rdata")

df1 <- data.frame( "B2" = eval_bqn2, "B10" = eval_bqn3, "L-B" = eval_lbqn, "SQ1" = eval_sq1, "ZAL" = eval_zal2, "DAAR" = eval_dar)

df2 <- data.frame("B2" = time_bqn2, "B10" = time_bqn3, "L-B" = time_lbqn, "SQ1" = time_sq1, "ZAL" = time_zal2, "DAAR" = time_dar)


pdf(file = "Out/quad-boxplot_sd1e3.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
boxplot(df1, xlab = "Acceleration Method", ylab = "No. of evaluations")
boxplot(df2, xlab = "Acceleration Method", ylab = "Time")
dev.off()
