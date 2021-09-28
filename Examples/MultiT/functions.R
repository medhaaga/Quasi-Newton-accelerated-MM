rm(list = ls())


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
    like = like - ((1+dim)/2) * log(1 + t(data[i,] - mu) %*% sig.inv %*% (data[i,] - mu))
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

param_constraint <- function(par)
{
  P <- length(par)
  p <- (-3 + sqrt(9 + (8*P)))/2
  vec <- par[-c(1,p)]
  mat <- VecToMat(vec, p)
  if(min(eigen(mat)$values) < 0)
    return(FALSE)
  else
    return(TRUE)
}
