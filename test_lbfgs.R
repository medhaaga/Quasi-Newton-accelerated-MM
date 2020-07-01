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
  sigma <- VecToMat(sigma, p)
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
  
  return (c(new.mu, upper.triangle(new.sigma, diag = TRUE)))
}

objective_fn <- function(current, p, P, data){
  new <- update(current[1:p], current[(p+1):P], data)
  objective <- new - current
  return (objective)
}

lbfgs_2loop <- function(t, m, m.u, m.v, current, p, P, data){
  
  next1 <- update(current[1:p], current[(p+1):P], data)
  next2 <- update(next1[1:p], next1[(p+1):P], data)
  
  u_t <- next1 - current
  v_t <- next2 - 2*next1 + current
  G.current <- u_t
  gamma_t <- dot(u_t, v_t)/dot(v_t, v_t)
  H_init <- gamma_t * diag(P)
  q <- G.current
  alpha <- rep(0,min(m, t-1))
  
  if (t >= 2){
    for (i in 1:min(m, (t-1))){
      rho <- 1/dot(m.v[[i]], m.v[[i]])
      alpha[i] <- rho*dot(m.v[[i]], q)
      q <- q - alpha[i]*m.v[[i]]
    }
  }
  
  r <- H_init %*% q
  if(t >= 2){
    for (i in 1:min(m, t-1)){
      r <- r + alpha[i]*m.u[[i]] 
    }
  }
  return (r)
}


##################################################

p <- 10
P <- 65
N <- 100
mu <- rep(0, p)
u <- matrix(rnorm(p*p), p, p)
sigma <- t(u) %*% u
data <- rmvc(n=N, mu = mu, S = sigma)


############################3
####L-BFGS for m=10
#############################

m <- 10
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
  start <- c(mu0, upper.triangle(sigma0, diag = TRUE))
  
  epsilon <- 1e-7
  current <- start
  theta <- start
  diff <- 100
  iter <- 0
  m.u <- list()
  m.v <- list()
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    
    direction <- -lbfgs_2loop(t = iter, m, m.u, m.v, current, p, P, data)
    foo <- current + direction
    
    theta <- list(foo[1:p], VecToMat(foo[(p+1):P], p))
    del.mu <- theta[[1]] - current[1:p]
    del.sigma <- theta[[2]] - VecToMat(current[(p+1):P], p)
    
    while(min(eigen(theta[[2]])$values) <= 0){
      theta[[2]] <- current[[2]] + del.sigma/2
      theta[[1]] <- current[[1]] + del.mu/2
      fails = fails+1
    }
    theta.vec <- c(theta[[1]], upper.triangle(theta[[2]], diag = TRUE))
    if(iter >=2){
      for (i in min(iter,m):2){
        m.u[[i]] <- m.u[[i-1]]
        m.v[[i]] <- m.v[[i-1]]
      }
    }
    m.u[[1]] <- objective_fn(theta.vec, p, P, data)
    m.v[[1]] <- objective_fn(m.u[[1]] + theta.vec, p, P, data)
    
    
    diff <- norm((theta.vec - current), type = "2")
    current <- theta.vec
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
