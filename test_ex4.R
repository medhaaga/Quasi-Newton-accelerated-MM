library(LaplacesDemon)
library(pracma)

VecToMat <- function(vec, p){
  a <- matrix(0, p, p)
  a[upper.tri(a, diag = TRUE)] <- vec
  a = a + t(a)
  diag(a) <- diag(a)/2
  return(a)
}

likelihood <- function(mu, sigma, data){
  
  n <- 100
  p <- 10
  like <- -(n*log(det(sigma)))/2
  sig.inv <- solve(sigma)
  for (i in 1:n){
    like = like - 5.5 * log(t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
  }
  return(like)
}

update <-  function(mu, sigma, data){
  n <- 100
  p <- 10
  sig.inv <- solve(sigma)
  weights <- as.matrix(rep(0, n))
  
  new.mu <- rep(0,p)
  new.sigma <- matrix(0, p, p)
  
  for (i in 1:n){
    weights[i] <- (1+p)/(1 + (t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu)))
    new.mu <- new.mu + weights[i]*as.matrix(data[i,])
    new.sigma <- new.sigma + weights[i]*((data[i,] - mu) %*% t(data[i,] - mu))
  }
  new.mu <- new.mu/sum(weights)
  new.sigma <- new.sigma/n
  
  return (list(new.mu, new.sigma))
}
##################################################

p <- 10``
P <- 65
N <- 100
mu <- rep(0, p)
u <- matrix(rnorm(p*p), p, p)
sigma <- t(u) %*% u
data <- rmvc(n=N, mu = mu, S = sigma)

###########################################
## Naive Algorithm
###########################################

n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0

for (j in 1:n){
  print(j)
  
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  theta <- start
  current <- theta
  diff <- 100
  epsilon <- 1e-7
  iter <- 0
  
  start.time <- Sys.time()
  while((diff > epsilon))
  {
    iter <- iter + 1
    if(iter %% 1000 == 0) print(iter)
    theta <- update(current[[1]], current[[2]], data)
    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))

##########################################
## Zhou's quasi-Newton for q=1
##########################################

n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-9
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0

  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    
    u <- as.matrix(next1 - c(current[[1]], upper.triangle(current[[2]], diag = TRUE)))
    v <- as.matrix(next2 - next1)
    foo1 <- -(dot(u,u))/(dot(u, (v-u)))
    foo2 <- (1 - foo1)*next1 + foo1*next2
    
    a <- matrix(0, p, p)
    a[upper.tri(a, diag = TRUE)] <- foo2[11:65]
    a = a + t(a)
    diag(a) <- diag(a)/2
    
    
    theta <- list(foo2[1:10], a)
    
    if(min(eigen(theta[[2]])$values) <= 0){
      print("Failed")
      fails = fails+1
      break
    }
    
    if(likelihood(theta[[1]], theta[[2]], data) < likelihood(current[[1]], current[[2]], data)){
      b <- matrix(0, p, p)
      b[upper.tri(b, diag = TRUE)] <- next2[11:65]
      b = b + t(b)
      diag(b) <- diag(b)/2
      theta <- list(next2[1:10], b)
      count = count+1
    }
    
    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)

########################################
## Classical BFGS
########################################


n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  
  H.current <- -diag(P)
  H.next <- H.current
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- (next1 - vec.current)
    u <- G.current
    v <- as.matrix(next2 - next1 - G.current)
    H.next <- (H.current %*% (diag(P) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next
    
    foo <- vec.current - H.current %*% G.current
    theta <- list(foo[1:10], VecToMat(foo[11:65], p))
    del.mu <- theta[[1]] - current[[1]]
    del.sigma <- theta[[2]] - current[[2]]
    
    while(min(eigen(theta[[2]])$values) <= 0){
      theta[[2]] <- current[[2]] + del.sigma/2
      theta[[1]] <- current[[1]] + del.mu/2
      fails = fails+1
    }

    
    #if(likelihood(theta[[1]], theta[[2]], data) < likelihood(current[[1]], current[[2]], data)){
    #  theta <- list(next2[1:10], VecToMat(next2[11:65], p))
    #  count = count+1
    #}
    
    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)


############################3
####BFGS + SQUAREM
#############################


n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  H.current <- -diag(P)
  H.next <- H.current
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- as.matrix(next1 - vec.current)
    u <- G.current
    v <- as.matrix(next2 - next1 - G.current)
    
    H.next <- (H.current %*% (diag(P) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next
    
    foo <- vec.current - H.current %*% G.current
    theta <- list(foo[1:10], VecToMat(foo[11:65], p))
    del.mu <- theta[[1]] - current[[1]]
    del.sigma <- theta[[2]] - current[[2]]
    
    while(min(eigen(theta[[2]])$values) <= 0){
      theta[[2]] <- current[[2]] + del.sigma/2
      theta[[1]] <- current[[1]] + del.mu/2
      fails = fails+1
    }

    #if(likelihood(theta[[1]], theta[[2]], data) < likelihood(current[[1]], current[[2]], data)){
    #  theta <- list(next2[1:10], VecToMat(next2[11:65], p))
    #  count = count+1
    #}
    
    H.current <- dot(u,v)/dot(v,v) * diag(P)
    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)


##########################################
#### SqS1
##########################################



n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- as.matrix(next1 - vec.current)
    u <- G.current
    v <- as.matrix(next2 - next1 - G.current)
    
    foo1 <- (dot(u,u))/(dot(u,v))
    foo2 <- vec.current - 2*foo1*u + (foo1^2)*v
    theta <- list(foo2[1:10], VecToMat(foo2[11:65], p))
    
    del.mu <- theta[[1]] - current[[1]]
    del.sigma <- theta[[2]] - current[[2]]
    
    while(min(eigen(theta[[2]])$values) <= 0){
      theta[[2]] <- current[[2]] + del.sigma/2
      theta[[1]] <- current[[1]] + del.mu/2
      fails = fails+1
    }
    
    diff <- diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)


##########################################
#### SqS2
##########################################


n <- 100
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- as.matrix(next1 - vec.current)
    u <- G.current
    v <- as.matrix(next2 - next1 - G.current)
    
    foo1 <- (dot(u,v))/(dot(v,v))
    foo2 <- vec.current - 2*foo1*u + (foo1^2)*v
    theta <- list(foo2[1:10], VecToMat(foo2[11:65], p))
    
    if(min(eigen(theta[[2]])$values) <= 0){
      print("Failed")
      fails = fails+1
      break
    }

    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)


##########################################
#### SqS3
##########################################



n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    next1 <- update(current[[1]], current[[2]], data)
    next2 <- update(next1[[1]], next1[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    next2 <- c(next2[[1]], upper.triangle(next2[[2]], diag = TRUE))
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- as.matrix(next1 - vec.current)
    u <- G.current
    v <- as.matrix(next2 - next1 - G.current)
    
    foo1 <- sqrt((dot(u,u))/(dot(v,v)))
    foo2 <- vec.current - 2*foo1*u + (foo1^2)*v
    theta <- list(foo2[1:10], VecToMat(foo2[11:65], p))
    
    if(min(eigen(theta[[2]])$values) <= 0){
      print("Failed")
      fails = fails+1
      break
    }

    diff <- diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
  }
  
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(count/sum(evals))
print(fails/n)

################################
#### JnJ-QN1
################################


n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- list()

epsilon <- 1e-7
fails <- 0

for (j in 1:n){
  print(j)
  data <- rmvc(n=N, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- list(mu0, sigma0)

  current <- start
  theta <- start
  
  H.current <- -diag(P)
  H.next <- H.current
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    next1 <- update(current[[1]], current[[2]], data)
    next1 <- c(next1[[1]], upper.triangle(next1[[2]], diag = TRUE))
    vec.current <- c(current[[1]], upper.triangle(current[[2]], diag = TRUE))
    G.current <- as.matrix(next1 - vec.current)
    
    foo <- vec.current - H.current %*% G.current
    theta <- list(foo[1:10], VecToMat(foo[11:65], p))
    del.sigma <- theta[[2]] - current[[2]]
    
    while(min(eigen(theta[[2]])$values) <= 0){
      theta[[2]] <- current[[2]] + del.sigma/2
      fails = fails+1
    }
    
    theta.next <- update(theta[[1]], theta[[2]], data)
    vec.theta <- c(theta[[1]], upper.triangle(theta[[2]], diag = TRUE))
    vec.theta.next <- c(theta.next[[1]], upper.triangle(theta.next[[2]], diag = TRUE))
    
    G.theta <- vec.theta.next - vec.theta
    del.theta <- vec.theta - vec.current
    del.G <- G.theta - G.current
    foo <- t(del.theta) %*% H.current
    H.next <- H.current + ((del.theta - (H.current %*% del.G)) %*% foo)/as.numeric(foo %*% del.G)
    H.current <- H.next
    
    
    diff <- norm(c(theta[[1]] - current[[1]], upper.triangle(theta[[2]] - current[[2]], diag = TRUE)), type = "2")
    current <- theta
    
  }
  end.time <- Sys.time()
  
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[[j]] <- theta
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(fails/n)
