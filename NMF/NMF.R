library(LaplacesDemon)
library(pracma)
library(SQUAREM)
library(BfgsQN)
source("qnamm.r")
set.seed(1)

frobenius <- function(now, V)
{
  dim <- round((length(now)/2)^.5)
  W.now <- matrix(now[1:dim^2], dim, dim)
  H.now <- matrix(now[-(1:dim^2)], dim, dim)
  temp <- (norm((V - (W.now %*% H.now)), type = "F"))
  return(temp)
}

KL_div <- function(now, V)
{
  r <- nrow(V)
  c <- ncol(V)
  dim <- round((length(now)/2)^.5)
  W.now <- matrix(now[1:dim^2], dim, dim)
  H.now <- matrix(now[-(1:dim^2)], dim, dim)

  A = V
  B = W.now %*% H.now
  sum <- 0
  for (i in 1:r)
  {
    for (j in 1:c)
    {
      sum <- sum + (A[i,j]*log(A[i,j]/B[i,j]) - A[i,j] + B[i,j])
    }
  }
  return(sum)
}

update1 <- function(now, V)
{
  dim <- (length(now)/2)^.5
  W.now <- matrix(now[1:dim^2], dim, dim)
  H.now <- matrix(now[-(1:dim^2)], dim, dim)

  H.new <- H.now * ((t(W.now) %*% V) / (t(W.now) %*% W.now %*% H.now))
  W.new <- W.now * ((V %*% t(H.new))/(W.now %*% H.new %*% t(H.new)))

  return (c(c(W.new), c(H.new)))
}


dim <- 10
V <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
A <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
B <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
start.all  <- matrix(rgamma(n= 2*dim*dim*N, shape = 5, scale = 1), nrow = N, ncol = 2*dim*dim)
tol <- 0.01

###############################################
#### MM Algorithm
###############################################

N <- 10
time_rep.mm <- rep(0,N)
fpevals.mm <- rep(0,N)
obj.values.mm <- rep(0, N)


for (j in 1:N)
{
  print(j)
  start  <- start.all[j,]

  now <- start
  new <- start
  iter <- 1
  diff <- 100


  start.time <- Sys.time()
  while(diff > tol)
  {
    if (iter %% 1000 == 0) print(diff)
    new <- update1(now, V)
    diff <- norm(new-now, type = "2")
    now <- new
    iter <- iter +1
  }
  end.time <- Sys.time()

  time_rep.mm [j] <- end.time - start.time
  fpevals.mm[j] <- iter
  obj.values.mm[j] <- frobenius(new, V)
}

print(quantile(time_rep.mm, probs = c(0, .5, 1)))
print(quantile(fpevals.mm, probs = c(0, .5, 1)))
print(quantile(obj.values.mm, probs = c(0, .5, 1)))

##################################################
#### BFGS
##################################################


N <- 10
time_rep.bfgs <- rep(0,N)
fevals.bfgs <- rep(0,N)
levals.bfgs <- rep(0, N)
obj.values.bfgs <- rep(0, N)
fails <- 0

for (j in 1:N)
{
  print(j)
  start  <- start.all[j,]

  start.time <- Sys.time()
  fp <- BFGS(par = start, V=V, fixptfn = update1, objfn = frobenius, control = list(tol = tol, objfn.inc = 1, maxiter = 5e4))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep.bfgs [j] <- end.time - start.time
    fevals.bfgs[j] <- fp$fpevals
    obj.values.bfgs[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.bfgs [j] <- NA
    fevals.bfgs[j] <- NA
    obj.values.bfgs[j] <- NA
  }
}

print(quantile(time_rep.bfgs, probs = c(0, .5, 1)))
print(quantile(fevals.bfgs, probs = c(0, .5, 1)))
print(quantile(obj.values.bfgs, probs = c(0, .5, 1)))



#####################################################
#### L-BFGS
#####################################################


N <- 10
time_rep.lbfgs <- rep(0,N)
fevals.lbfgs <- rep(0,N)
levals.lbfgs <- rep(0, N)
obj.values.lbfgs <- rep(0, N)
tol = 1e-3
fails <- 0

for (j in 1:N)
{
  print(j)
  start  <- start.all[j,]


  start.time <- Sys.time()
  fp <- LBFGS(par = start, V=V, fixptfn = update1, objfn = frobenius, control = list(tol = tol, objfn.inc = 1, maxiter = 5e4))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep.lbfgs [j] <- end.time - start.time
    fevals.lbfgs[j] <- fp$fpevals
    levals.lbfgs[j] <- fp$objfevals
    obj.values.lbfgs[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.lbfgs [j] <- NA
    fevals.lbfgs[j] <- NA
    levals.lbfgs[j] <- NA
    obj.values.lbfgs[j] <- NA
  }
}

print(quantile(time_rep.lbfgs, probs = c(0, .5, 1)))
print(quantile(fevals.lbfgs, probs = c(0, .5, 1)))
print(quantile(obj.values.lbfgs, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 1
##########################################

N <- 10
time_rep.sq1 <- rep(0,N)
fevals.sq1 <- rep(0,N)
levals.sq1 <- rep(0, N)
obj.values.sq1 <- rep(0, N)
tol = 1e-3
fails <- 0


for (j in 1:N){
  print(j)
  start  <- start.all[j,]

  start.time <- Sys.time()
  fp <- squarem(start, V=V, fixptfn = update1, objfn = frobenius, control = list(K=1, tol = tol, method = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fevals.sq1[j] <- fp$fpevals
    levals.sq1[j] <- fp$objfevals
    obj.values.sq1[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep.sq1 [j] <- NA
    fevals.sq1[j] <- NA
    levals.sq1[j] <- NA
    obj.values.sq1[j] <- NA
  }
}
print(quantile(time_rep.sq1, probs = c(0, .5, 1)))
print(quantile(fevals.sq1, probs = c(0, .5, 1)))
print(quantile(levals.sq1, probs = c(0, .5, 1)))
print(quantile(obj.values.sq1, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 2
##########################################

N <- 10
time_rep.sq2 <- rep(0,N)
fevals.sq2 <- rep(0,N)
levals.sq2 <- rep(0, N)
obj.values.sq2 <- rep(0, N)
tol = 1e-3
fails <- 0


for (j in 1:N){
  print(j)
  start  <- start.all[j,]

  start.time <- Sys.time()
  fp <- squarem(start, V=V, fixptfn = update1, objfn = frobenius, control = list(K=1, tol = tol, method = 2))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fevals[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.sq2[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    fevals[j] <- NA
    levals[j] <- NA
    obj.values.sq2[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fevals, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.sq2, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 3
##########################################

N <- 10
time_rep <- rep(0,N)
fevals <- rep(0,N)
levals <- rep(0, N)
obj.values.sq3 <- rep(0, N)
tol = 1e-3
fails <- 0


for (j in 1:N){
  print(j)
  A <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
  B <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
  start  <- c(c(A), c(B))

  start.time <- Sys.time()
  fp <- squarem(start, V=V, fixptfn = update1, objfn = frobenius, control = list(K=1, tol = tol, method = 3))

  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fevals[j] <- fp$fpevals
    levals[j] <- fp$objfevals
    obj.values.sq3[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    fevals[j] <- NA
    levals[j] <- NA
    obj.values.sq3[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fevals, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.sq3, probs = c(0, .5, 1)))

#################################################
##### ZAL; q=2
#################################################

q <- 2
N <- 10
time_rep <- rep(0,N)
fevals <- rep(0,N)
levals <- rep(0,N)
obj.values.zal <- rep(0, N)
accept <- rep(0, N)
reject <- rep(0, N)
tol = 1e-3
fails <- 0


for (j in 1:N){
  print(j)
  A <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
  B <- matrix(rgamma(n= dim*dim, shape = 5, scale = 1), dim, dim)
  start  <- c(c(A), c(B))

  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update1, qn=q, fx_obj = frobenius, max_iter = 5e4, tol = tol, V=V)
  end.time <- Sys.time()

  if (fp$convergence == TRUE){
  time_rep[j] <- end.time - start.time
  fevals[j] <- fp$fevals
  levals[j] <- fp$levals
  obj.values.zal[j] <- fp$objective
  accept[j] <- fp$accept
  reject[j] <- fp$reject}
  else{
    time_rep[j] <- NA
    fevals[j] <- NA
    levals[j] <- NA
    obj.values.zal[j] <- NA
    accept[j] <- NA
    reject[j] <- NA
  }
}

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fevals, probs = c(0, .5, 1)))
print(quantile(levals, probs = c(0, .5, 1)))
print(quantile(obj.values.zal, probs = c(0, .5, 1)))
print(quantile(accept, probs = c(0, .5, 1)))
print(quantile(reject, probs = c(0, .5, 1)))


####################################
#### Boxplots
#####################################

####################################
#### Boxplots and scatterplots
#####################################



df1 <- data.frame("MM" = obj.values.mm, "SqS1" = obj.values.sq1, "SqS2" = obj.values.sq2,
                  "SqS3" = obj.values.sq3, "BFGS" = obj.values.bfgs, "LBFGS" = obj.values.lbfgs)

df2 <- data.frame("MM" = fevals.mm, "SqS1" = fevals.sq1, "SqS2" = fevals.sq2,
                  "SqS3" = fevals.sq3, "BFGS" = fevals.bfgs, "LBFGS" = fevals.lbfgs)

df3 <- data.frame("MM" = time_rep.mm, "SqS1" = time_rep.sq1, "SqS2" = time_rep.sq2,
                  "SqS3" = time_rep.sq3, "BFGS" = time_rep.bfgs, "LBFGS" = time_rep.lbfgs)


pdf(file = "Out/boxplot_NMF.pdf", width = 12, height = 6)
par(mfrow = c(1,3))
boxplot(df1, xlab = "Algorithm", ylab = "Objective value")
boxplot(df2, xlab = "Algorithm", ylab = "No. of evaluations")
boxplot(df3, xlab = "Algorithm", ylab = "Time")
dev.off()

pdf(file = "Out/scatter_plot_NMF.pdf")
plot(fevals.mm, obj.values.mm, pch = 19, col = "cadetblue1", xlab = "Number of updates", ylab = "Objective Value",
     xlim = range(fevals.mm, fevals.bfgs, fevals.lbfgs,fevals.sq1, fevals.sq2, fevals.sq3),
     ylim = range(obj.values.bfgs, obj.values.lbfgs, obj.values.mm, obj.values.sq1, obj.values.sq2, obj.values.sq3))
points(fevals.bfgs, obj.values.bfgs, pch = 19, col = "cyan3")
points(fevals.lbfgs, obj.values.lbfgs, pch = 19, col = "hotpink")
points(fevals.sq1, obj.values.sq1, pch = 19, col = "lightpink")
points(fevals.sq2, obj.values.sq2, pch = 19, col = "coral")
points(fevals.sq3, obj.values.sq3, pch = 19, col = "blue")
legend("bottomright", legend = c("MM", "BFGS", "L-BFGS", "SqS1", "SqS2", "SqS3"), col = c("cyan3", "hotpink", "lightpink", "coral", "blue"), pch = 19)
dev.off()
