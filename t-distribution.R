library(LaplacesDemon)
library(pracma)
library(SQUAREM)
library(BfgsQN)

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
  like <- -(n*log(det(sigma)))/2
  sig.inv <- solve(sigma)
  for (i in 1:n){
    like = like - ((1+dim)/2) * log(t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
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

dim <- 50

P <- (dim/2)*(dim+3)
n <- 100
mu <- rep(0, dim)
u <- matrix(rnorm(dim*dim), dim, dim)
sigma <- t(u) %*% u

###########################################
## EM Algorithm
###########################################

N <- 1
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0

for (j in 1:N){
  print(j)

  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag=  TRUE))

  start.time <- Sys.time()
  fp <- fpiter(par = start, n=n, dim=dim, data=data, fixptfn = update, objfn = likelihood, control = list(tol = epsilon, maxiter = 1e5))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))

###########################################
## PX-EM
###########################################

N <- 1
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0

for (j in 1:N){
  print(j)

  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag=  TRUE))

  start.time <- Sys.time()
  fp <- fpiter(par = start, n=n, dim=dim, data=data, fixptfn = update_pxem, objfn = likelihood, control = list(tol = epsilon, maxiter = 1e3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))

##########################################
## ZAL, q=1
##########################################

N <- 1
time_rep <- rep(0,N)
fpevals <- rep(0,N)
obj.values <- rep(0, N)
fallback.prop <- rep(0, N)


for (j in 1:N){
  print(j)
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag = TRUE))

  epsilon <- 1e-7
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0
  nrep <- 0
  fallback <- 0

  start.time <- Sys.time()
  while(diff > epsilon){

    nrep = nrep + 1
    next1 <- update(current, n, dim, data)
    next2 <- update(next1, n, dim, data)
    iter <- iter + 2

    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1)
    alpha <- -as.numeric(crossprod(u))/as.numeric(crossprod(u, (v-u)))
    theta <- (1 - alpha)*next1 + alpha*next2

    sig.old <- VecToMat(current[-(1:dim)], dim)
    sig.new <- VecToMat(theta[-(1:dim)], dim)

    del.mu <- theta[1:dim] - current[1:dim]
    del.sigma <- sig.new - sig.old

    while(min(eigen(sig.new)$values) <= 0){
      sig.new <- sig.old + del.sigma/2
      theta[1:dim] <- current[1:dim] + del.mu/2
      del.sigma <- del.sigma/2
      del.mu <- del.mu/2
      print("fail")
    }

    theta[-(1:dim)] <- upper.triangle(sig.new, dim)


    if(likelihood(theta, n, dim, data) < likelihood(current, n, dim, data)){
      theta <- next2
      fallback = fallback+1
    }

    diff <- norm(theta - current, type = "2")
    current <- theta
  }
  end.time <- Sys.time()

  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  obj.values[j] <- likelihood(theta, n, dim, data)
  fallback.prop[j] <- fallback/nrep

}

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(quantile(fallback.prop, probs = c(0, .5, 1)))


########################################
## Classical BFGS
########################################

N <- 1
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0
epsilon <- 1e-7

for (j in 1:N){
  print(paste("nrep =", j))
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag = TRUE))

  start.time <- Sys.time()
  fp <- BFGS(par = start, dim = dim, n= n, data = data, fixptfn = update, objfn = likelihood, control = list(tol = epsilon, objfn.inc = 20))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }

}

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))


##########################################
#### SQUAREM - 1
##########################################

N <- 100
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0

for (j in 1:N){
  print(j)

  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag=  TRUE))

  start.time <- Sys.time()
  fp <- squarem(start, dim = dim, n= n, data = data, fixptfn = update, objfn = likelihood, control = list(K=1, tol = epsilon, method = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 2
##########################################

N <- 100
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0

for (j in 1:N){
  print(j)

  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag=  TRUE))

  start.time <- Sys.time()
  fp <- squarem(start, dim = dim, n= n, data = data, fixptfn = update, objfn = likelihood, control = list(K=1, tol = epsilon, method = 2))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))

##########################################
#### SQUAREM - 3
##########################################

N <- 100
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)
fails <- 0

for (j in 1:N){
  print(j)

  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag=  TRUE))

  start.time <- Sys.time()
  fp <- squarem(start, dim = dim, n= n, data = data, fixptfn = update, objfn = likelihood, control = list(K=1, tol = epsilon, method = 3))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    evals[j] <- fp$fpevals
    obj.values[j] <- fp$value.objfn
  } else{
    fails = fails+1
    time_rep [j] <- NA
    evals[j] <- NA
    obj.values[j] <- NA
  }
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))

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
