library(LaplacesDemon)
library(pracma)

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
    like = like - 5.5 * log(t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
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


objective_fn <- function(current, n, dim, data){
  new <- update(current, n, dim, data)
  objective <- new - current
  return (objective)
}


lbfgs_2loop <- function(t, m, m.u, m.v, current, n, dim, data){
  
  P <- (dim/2)*(dim+3)
  next1 <- update(current, n, dim, data)
  next2 <- update(next1, n, dim, data)
  
  u_t <- next1 - current
  v_t <- next2 - 2*next1 + current
  G.current <- u_t
  gamma_t <- as.numeric(crossprod(u_t, v_t))/as.numeric(crossprod(v_t, v_t))
  H_init <- gamma_t*diag(P)
  q <- G.current
  alpha <- rep(0,min(m, t-1))
  
  if (t >= 2){
    for (i in 1:min(m, (t-1))){
      rho <- 1/as.numeric(crossprod(m.v[[i]], m.v[[i]]))
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

dim <- 50
P <- (dim/2)*(dim+3)
n <- 1000
mu <- rep(0, dim)
u <- matrix(rnorm(dim*dim), dim, dim)
sigma <- t(u) %*% u


############################3
####L-BFGS for m=10
#############################

m <- 10
N <- 1
epsilon <- 1e-7
time_rep <- rep(0,N)
evals <- rep(0,N)
obj.values <- rep(0, N)


for (j in 1:N){
  print(j)
  data <- rmvc(n=n, mu = mu, S = sigma)
  mu0 <- colMeans(data)
  sigma0 <- cov(data)
  start <- c(mu0, upper.triangle(sigma0, diag = TRUE))
  
  current <- start
  theta <- start
  diff <- 100
  iter <- 0
  levals <- 0
  m.u <- list()
  m.v <- list()
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    levals <- levals+1
    direction <- -lbfgs_2loop(t = levals, m, m.u, m.v, current, n, dim, data)
    iter <- iter+2
    theta <- current + direction
    
    sig.new <- VecToMat(theta[-(1:dim)], dim)
    sig.old <- VecToMat(current[-(1:dim)], dim)
    del.mu <- theta[1:dim] - current[1:dim]
    del.sigma <-  sig.new - sig.old
    
    while(min(eigen(sig.new)$values) <= 0){
      sig.new <- sig.old + del.sigma/2
      theta[1:dim] <- current[1:dim] + del.mu/2
      del.sigma <- del.sigma/2
      del.mu <- del.mu/2
      print("fail")
    }
    
    theta[-(1:dim)] <- upper.triangle(sig.new, diag = TRUE)
    
    if(levals >=2){
      for (i in min(levals,m):2){
        m.u[[i]] <- m.u[[i-1]]
        m.v[[i]] <- m.v[[i-1]]
      }
    }
    m.u[[1]] <- objective_fn(theta, n, dim, data)
    m.v[[1]] <- objective_fn(m.u[[1]] + theta, n, dim, data) - m.u[[1]]
    
    
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
