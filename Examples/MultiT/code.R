library(LaplacesDemon)
library(pracma)
library(SQUAREM)
library(BfgsQN)
source("qnamm.r")

VecToMat <- function(vec, dim){
  a <- matrix(0, dim, dim)
  a[upper.tri(a, diag = TRUE)] <- vec
  a = a + t(a)
  diag(a) <- diag(a)/2
  return(a)
}

likelihood <- function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  like <- (n*log(det(sigma)))/2
  sig.inv <- solve(sigma)
  for (i in 1:n){
    like = like + ((1+dim)/2) * log(1 + t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
  }
  return(like)
}

update <-  function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  sig.inv <- solve(sigma)
  weights <- as.matrix(rep(0, n))

  new.mu <- rep(0,dim)
  new.sigma <- matrix(0, dim, dim)

  for (i in 1:n){
    weights[i] <- (1+dim)/(1 + (t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu)))
    new.mu <- new.mu + weights[i]*as.matrix(data[i,])
    new.sigma <- new.sigma + weights[i]*((data[i,] - mu) %*% t(data[i,] - mu))
  }
  new.mu <- new.mu/sum(weights)
  new.sigma <- new.sigma/n

  return (c(new.mu, upper.triangle(new.sigma, diag = TRUE)))
}

update_pxem <-  function(par, n, dim, data){
  mu <- par[1:dim]
  sigma <- VecToMat(par[-(1:dim)], dim)
  sig.inv <- solve(sigma)
  weights <- as.matrix(rep(0, n))

  new.mu <- rep(0,dim)
  new.sigma <- matrix(0, dim, dim)

  for (i in 1:n){
    weights[i] <- (1+dim)/(1 + (t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu)))
    new.mu <- new.mu + weights[i]*as.matrix(data[i,])
    new.sigma <- new.sigma + weights[i]*((data[i,] - mu) %*% t(data[i,] - mu))
  }
  new.mu <- new.mu/sum(weights)
  new.sigma <- new.sigma/sum(weights)

  return (c(new.mu, upper.triangle(new.sigma, diag = TRUE)))
}


##################################################

set.seed(10)
dim <- 10
tol <- 1e-5
P <- (dim/2)*(dim+3)
n <- 100
mu <- rep(0, dim)
u <- matrix(rnorm(dim*dim), dim, dim)
sigma <- t(u) %*% u
N <- 1
start_rep <- matrix(0, nrow = N, ncol = P)
for (i in 1:N)
{
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start_rep[i,] <- c(mu0, upper.triangle(sigma0, diag=  TRUE))
}


###########################################
## EM Algorithm
###########################################

time_mm <- rep(0, N)
obj_mm <- rep(0, N)
eval_mm <- rep(0, N)

for (i in 1:N){
  print(i)
  now <- start_rep[i,]
  new <- start_rep[i,]
  objval <- likelihood(now, n=n, dim=dim, data=data)
  chain_mm <- objval
  diff <- 100
  iter <- 0
  start.time <- Sys.time()
  chain <- matrix(0, nrow = 1e6, ncol = P)
  while((diff > tol))
  {
    iter <- iter + 1
    if(iter %% 1000 == 0) print(iter)
    new <- update(now, n=n, dim=dim, data=data)
    objval <- likelihood(new, n=n, dim=dim, data=data)
    chain_mm <- c(chain_mm, objval)
    chain[iter,] <- new
    diff <- sqrt(crossprod(new-now))
    now <- new
  }
  end.time <- Sys.time()
  chain <- chain[1:iter,]
  time_mm[i] <- end.time - start.time
  obj_mm[i] <- likelihood(new, n=n, dim=dim, data=data)
  eval_mm[i] <- iter

  #plot(chain[,159], chain[,160], xlim = range(chain[,159]), ylim = range(chain[,160]))
  #points(chain[1,159], chain[1,160], col = "green", pch=19, cex=1.5)
  #points(chain[iter,159], chain[iter,160], col = "red", pch=19, cex=1.5)

}

print(quantile(time_mm, probs = c(0, .5, 1)))
print(quantile(eval_mm, probs = c(0, .5, 1)))
print(quantile(obj_mm, probs = c(0, .5, 1)))

###########################################
## PX-EM
###########################################

time_pxem <- rep(0, N)
obj_pxem <- rep(0, N)
eval_pxem <- rep(0, N)

for (i in 1:N){
  print(i)
  start <- start_rep[i,]

  start.time <- Sys.time()
  fp <- fpiter(par = start, n=n, dim=dim, data=data, fixptfn = update_pxem,
               objfn = likelihood, control = list(tol = tol, maxiter = 1e5))
  end.time <- Sys.time()

  if(fp$convergence){
    time_pxem[i] <- end.time - start.time
    obj_pxem[i] <- fp$value.objfn
    eval_pxem[i] <- fp$fpevals
  } else{
    time_pxem[i] <- NA
    obj_pxem[i] <- NA
    eval_pxem[i] <- NA
  }
}

print(quantile(time_pxem, probs = c(0, .5, 1)))
print(quantile(eval_pxem, probs = c(0, .5, 1)))
print(quantile(obj_pxem, probs = c(0, .5, 1)))

##########################################
## ZAL, q=1
##########################################

time_zal <- rep(NA, N)
obj_zal <- rep(NA, N)
eval_zal <- rep(NA, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- qnamm(x = start, fx_mm = update, fx_obj = likelihood, qn=2, n=n, dim=dim, data=data, max_iter = 1e4, tol=tol)
  end.time <- Sys.time()

  time_zal[i] <- end.time - start.time
  obj_zal[i] <- fp$objective
  eval_zal[i] <- fp$fevals
  chain_zal <- fp$obj_iter
}

fails <- sum(is.na(obj_zal))
time_zal <- time_zal[!is.na(time_zal)]
eval_zal <- eval_zal[!is.na(eval_zal)]
obj_zal <- obj_zal[!is.na(obj_zal)]
print(quantile(time_zal, probs = c(0, .5, 1)))
print(quantile(eval_zal, probs = c(0, .5, 1)))
print(quantile(obj_zal, probs = c(0, .5, 1)))

########################################
## Classical BFGS1
########################################

time_bfgs1 <- rep(0, N)
obj_bfgs1 <- rep(0, N)
eval_bfgs1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BFGS(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
             control = list(qn=1, tol = tol, objfn.inc = 1, step.max = 1e7, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_bfgs1[i] <- end.time - start.time
  obj_bfgs1[i] <- fp$value.objfn
  eval_bfgs1[i] <- fp$fpevals
  chain_bfgs1 <- fp$p.inter[,(P+1)]
}

print(quantile(time_bfgs1, probs = c(0, .5, 1)))
print(quantile(eval_bfgs1, probs = c(0, .5, 1)))
print(quantile(obj_bfgs1, probs = c(0, .5, 1)))

########################################
## Classical bfgs2
########################################

time_bfgs2 <- rep(0, N)
obj_bfgs2 <- rep(0, N)
eval_bfgs2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- BFGS(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
             control = list(qn=2, tol = tol, objfn.inc = 1, step.max = 1e6, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  
  time_bfgs2[i] <- end.time - start.time
  obj_bfgs2[i] <- fp$value.objfn
  eval_bfgs2[i] <- fp$fpevals
  chain_bfgs2 <- fp$p.inter[,(P+1)]
}

print(quantile(time_bfgs2, probs = c(0, .5, 1)))
print(quantile(eval_bfgs2, probs = c(0, .5, 1)))
print(quantile(obj_bfgs2, probs = c(0, .5, 1)))

########################################
## L-BFGS
########################################

time_lbfgs <- rep(0, N)
obj_lbfgs <- rep(0, N)
eval_lbfgs <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- LBFGS(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
              control = list(m=10, tol = tol, objfn.inc = 10, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()
  time_lbfgs[i] <- end.time - start.time
  obj_lbfgs[i] <- fp$value.objfn
  eval_lbfgs[i] <- fp$fpevals
  chain_lbfgs <- fp$p.inter[,(P+1)]
  }
print(quantile(time_lbfgs, probs = c(0, .5, 1)))
print(quantile(eval_lbfgs, probs = c(0, .5, 1)))
print(quantile(obj_lbfgs, probs = c(0, .5, 1)))


  ##########################################
### SqS1
#############################################

time_sq1 <- rep(0, N)
obj_sq1 <- rep(0, N)
eval_sq1 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 1, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq1[i] <- end.time - start.time
  obj_sq1[i] <- fp$value.objfn
  eval_sq1[i] <- fp$fpevals
  chain_sq1 <- fp$p.inter[,(P+1)]
}

print(quantile(time_sq1, probs = c(0, .5, 1)))
print(quantile(eval_sq1, probs = c(0, .5, 1)))
print(quantile(obj_sq1, probs = c(0, .5, 1)))

##########################################
### SqS2
#############################################

time_sq2 <- rep(0, N)
obj_sq2 <- rep(0, N)
eval_sq2 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood, n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 2, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq2[i] <- end.time - start.time
  obj_sq2[i] <- fp$value.objfn
  eval_sq2[i] <- fp$fpevals
  chain_sq2 <- fp$p.inter[,(P+1)]
}
print(quantile(time_sq2, probs = c(0, .5, 1)))
print(quantile(eval_sq2, probs = c(0, .5, 1)))
print(quantile(obj_sq2, probs = c(0, .5, 1)))

##########################################
### SqS3
#############################################

time_sq3 <- rep(0, N)
obj_sq3 <- rep(0, N)
eval_sq3 <- rep(0, N)

for (i in 1:N){
  start <- start_rep[i,]
  start.time <- Sys.time()
  fp <- squarem(par = start, fixptfn = update, objfn = likelihood,n=n, dim=dim, data=data,
                control = list(K=1, tol = tol, method = 3, maxiter = 5e4, intermed = TRUE))
  end.time <- Sys.time()

  time_sq3[i] <- end.time - start.time
  obj_sq3[i] <- fp$value.objfn
  eval_sq3[i] <- fp$fpevals
  chain_sq3 <- fp$p.inter[,(P+1)]
}

print(quantile(time_sq3, probs = c(0, .5, 1)))
print(quantile(eval_sq3, probs = c(0, .5, 1)))
print(quantile(obj_sq3, probs = c(0, .5, 1)))

################################
#### JnJ-QN1
################################


N <- 100
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)

epsilon <- 1e-7

for (j in 1:N){

  print(j)
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag = TRUE))

  current <- start
  theta <- start

  H.current <- -diag(P)
  H.next <- H.current
  diff <- 100
  iter <- 0
  nrep <- 0
  start.time <- Sys.time()

  while(diff > epsilon){

    next1 <- update(current, n, dim, data)
    iter <- iter + 1

    G.current <- as.matrix(next1 - current)

    theta <- current - H.current %*% G.current

    sig.new <- VecToMat(theta[-(1:dim)], dim)
    sig.old <- VecToMat(current[-(1:dim)], dim)
    del.mu <- theta[1:dim] - current[1:dim]
    del.sigma <-  sig.new - sig.old

    while(min(eigen(sig.new)$values) <= 0){
      sig.new <- sig.old + del.sigma/2
      theta[1:dim] <- current[1:dim] + del.mu/2
      del.sigma <- del.sigma/2
      del.mu <- del.mu/2
    }

    theta[-(1:dim)] <- upper.triangle(sig.new, diag = TRUE)

    theta.next <- update(theta, n, dim, data)
    iter <- iter+1

    G.theta <- theta.next - theta
    del.theta <- theta - current
    del.G <- G.theta - G.current
    foo <- t(del.theta) %*% H.current
    H.next <- H.current + ((del.theta - (H.current %*% del.G)) %*% foo)/as.numeric(foo %*% del.G)
    H.current <- H.next


    diff <- norm((theta - current), type = "2")
    current <- theta

  }
  end.time <- Sys.time()

  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  obj.values[j] <- likelihood(theta, n, dim, data)

}

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))


save(time_mm, time_bfgs1, time_bfgs2, time_lbfgs, time_sq1, time_sq2, time_sq3, time_zal,
     eval_mm, eval_bfgs1, eval_bfgs2, eval_lbfgs, eval_sq1, eval_sq2, eval_sq3, eval_zal,
     obj_mm, obj_bfgs1, obj_bfgs2, obj_lbfgs, obj_sq1, obj_sq2, obj_sq3, obj_zal,
     c, chain_bfgs1, chain_bfgs2, chain_lbfgs, chain_sq1, chain_sq2, chain_sq3, chain_zal, file = "Out/multiT-objects.Rdata")

load(file = "Out/multiT-objects.Rdata")
pdf(file = "Out/multiT-running.pdf")
plot(1:50, chain_bfgs2[1:50], type = 'l', col = "lightskyblue", lwd=2, 
     xlab = "Ïterations", ylab = "Negative likelihood")
lines(1:50, chain_bfgs1[1:50], col = "red", lty=2, lwd=2)
lines(1:50, chain_lbfgs[1:50], col = "dodgerblue4", lty=3, lwd=2)
lines(1:50, chain_zal[1:50], col = "green", lty=4, lwd=2)
lines(1:50, chain_sq1[1:50], col = "pink", lty=5, lwd=2)
legend("topright", legend = c("BFGS1", "BFGS2", "LBFGS", "ZAL", "SQ1"), lty = 1:5, col = c("red", "lightskyblue", "dodgerblue4", "green", "pink"), lwd=2)
dev.off()
