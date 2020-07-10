set.seed(100)
library(pracma)


factorial <- function(num){
  fac = 1
  if(num < 0) {
    print("Sorry, factorial does not exist for negative numbers")
  } else if(num == 0) {
    return(1)
  } else {
    for(i in 1:num) {
      fac = fac * i
    }
    return(fac)
  }
}


likelihood <- function(data, p, mu1, mu2)
{
  like <- 0
  for (i in 0:9){
    like <- like + data[i+1,2] * log(((p * exp(-mu1) * (mu1 ^ i))/factorial(i)) + (((1-p) * exp(-mu2) * (mu2 ^ i))/factorial(i)))
  }
  return(like)
}

pi_hat <- function(p, mu1, mu2){
  
  pi <- matrix(0, nrow = 10, ncol = 2)
  denom <- rep(0, 10)
  for (i in 0:9){
    denom[i+1] <- (p * (mu1 ^ i) * exp(-mu1)) + ((1-p) * (mu2 ^ i) * exp(-mu2))
    pi[i+1,1] <- (p * (mu1 ^ i) * exp(-mu1))/denom[i+1]
    pi[i+1,2] <- 1 - pi[i+1,1]
  }
  return (pi)
}

update <- function(p, mu1, mu2, pi_hat, data){
  new_p <- dot(data[,2], pi_hat[,1])/sum(data[,2])
  new_mu1 <- dot((data[,1]*data[,2]), pi_hat[,1])/dot(data[,2], pi_hat[,1])
  new_mu2 <- dot((data[,1]*data[,2]), pi_hat[,2])/dot(data[,2], pi_hat[,2])
  return(c(new_p, new_mu1, new_mu2))
}

objective_fn <- function(current, p, data){
  
  pi <- pi_hat(current[1], current[2], current[3])
  next1 <- update(current[1], current[2], current[3], pi, data)
  objective <- next1 - current
  return (objective)
}

lbfgs_2loop <- function(t, m, m.u, m.v, current, p, data){
  
  pi <- pi_hat(current[1], current[2], current[3])
  next1 <- update(current[1], current[2], current[3], pi, data)
  pi <- pi_hat(next1[1], next1[2], next1[3])
  next2 <- update(next1[1], next1[2], next1[3], pi, data)
  
  u_t <- next1 - current
  v_t <- next2 - 2*next1 + current
  G.current <- u_t
  gamma_t <- dot(u_t, v_t)/dot(v_t, v_t)
  H_init <- gamma_t*diag(p)
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

############################################################3

deaths <- (0:9)
freq <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
data <- as.matrix(cbind(deaths, freq))
truth <- c(.3599, 1.256, 2.663)
#start <- c(.3, 1, 2.5)
p <- 3


############################3
####L-BFGS for m=10
#############################

m <- 10
fails <- 0

  start <- c(0.3, 1, 2.5)
  
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
    
    direction <- -lbfgs_2loop(t = iter, m, m.u, m.v, current, p, data)
    theta <- current + direction
    
    del <- theta - current
    
    while(theta[1]<=0.05 || theta[1]>= 0.95){
      theta <- current + del/2
      del <- del/2
      fails = fails+1
    }
    
    if(iter >=2){
      for (i in min(iter,m):2){
        m.u[[i]] <- m.u[[i-1]]
        m.v[[i]] <- m.v[[i-1]]
      }
    }
    m.u[[1]] <- objective_fn(theta, p, data)
    m.v[[1]] <- objective_fn(m.u[[1]] + theta, p, data) - m.u[[1]]
    
    
    diff <- norm((theta - current), type = "2")
    current <- theta
  }
  end.time <- Sys.time()
  
print(end.time - start.time)
print(iter)
print(theta)
  


