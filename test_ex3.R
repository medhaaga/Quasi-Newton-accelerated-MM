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

############################################################3

deaths <- (0:9)
freq <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
data <- as.matrix(cbind(deaths, freq))
truth <- c(.3599, 1.256, 2.663)
#start <- c(.3, 1, 2.5)
p <- 3

###########################################
## Naive Algorithm
###########################################

n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0
for (j in 1:n){

start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))

theta <- start
current <- theta
diff <- 100
epsilon <- 1e-9
iter <- 0

start.time <- Sys.time()
while((diff > epsilon))
{
  iter <- iter + 1
  if(iter %% 1000 == 0) print(iter)
  pi <- pi_hat(current[1], current[2], current[3])
  theta <- update(current[1], current[2], current[3], pi, data)
  diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
              likelihood(data, current[1], current[2], current[3]))/
          (abs(likelihood(data, current[1], current[2], current[3])) + 1)
  current <- theta
}
end.time <- Sys.time()
if (theta[1]< 0.05 || theta[1] > 0.95){
  fails = fails+1
  time_rep[j] <- NA
  evals[j] <- NA
  estimates[j,] <- NA
}
else{
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[j,] <- theta
}
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)


##########################################
## Zhou's quasi-Newton for q=1
##########################################

n <- 100
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))

  epsilon <- 1e-9
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0
  count <- 0

  start.time <- Sys.time()
  while(diff > epsilon){

    iter <- iter + 1
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1)
    foo <- -(dot(u,u))/(dot(u, (v-u)))
    theta <- (1 - foo)*next1 + foo*next2
    #if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
    #  print("Failed")
    #  theta <- NA
    #  iter <- NA
    #  break
    #}

    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
      count = count+1
    }
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    current <- theta
}
end.time <- Sys.time()

if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
  fails = fails+1
  time_rep[j] <- NA
  evals[j] <- NA
  estimates[j,] <- NA
}
else{
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[j,] <- theta
}

}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)


########################################
## Classical BFGS
########################################


n <- 100
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))

  epsilon <- 1e-9
  current <- start
  theta <- start
  
  H.current <- -diag(p)
  H.next <- H.current
  diff <- 100
  iter <- 0

  start.time <- Sys.time()
  while(diff > epsilon){

    iter <- iter + 1
    print(iter)
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)

    G.current <- next1 - current
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1 - G.current)

    H.next <- (H.current %*% (diag(p) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next


    theta <- current - H.current %*% G.current


    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      theta <- NA
      iter <- NA
      break
    }
    
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
      count = count+1
    }

    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    
    current <- theta
  }
  end.time <- Sys.time()
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  }
  else{
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)

##########################################
#### SqS1
##########################################



n <- 100
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
  
  epsilon <- 1e-9
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0
  count <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
  
    iter <- iter + 1
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)
  
    G.current <- next1 - current
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1 - G.current)
  
    foo <- (dot(u,u))/(dot(u,v))
  
    theta <- current - 2*foo*u + (foo^2)*v
    
    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      theta <- NA
      iter <- NA
      break
    }
    
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
      count = count+1
    }
    
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    current <- theta
  }
  
  end.time <- Sys.time()
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  }
  else{
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)


##########################################
#### SqS2
##########################################



n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 100), 2), round(runif(1, 1, 100)))
  
  epsilon <- 1e-9
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)
    
    G.current <- next1 - current
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1 - G.current)
    
    foo <- (dot(u,v))/(dot(v,v))
    
    theta <- current - 2*foo*u + (foo^2)*v
    
    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      theta <- NA
      iter <- NA
      break
    }
    
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
    }
    
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    current <- theta
  }
  
  end.time <- Sys.time()
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  }
  else{
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)


##########################################
#### SqS3
##########################################

n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 100), 2), round(runif(1, 1, 100)))
  
  epsilon <- 1e-9
  theta <- start
  current <- theta
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  
  while(diff > epsilon){
    
    iter <- iter + 1
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)
    
    G.current <- next1 - current
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1 - G.current)
    
    foo <- sqrt((dot(u,u))/(dot(v,v)))
    
    theta <- current - 2*foo*u + (foo^2)*v
    
    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      theta <- NA
      iter <- NA
      break
    }
    
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
    }
    
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    current <- theta
  }
  
  end.time <- Sys.time()
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  }
  else{
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)

############################3
####BFGS + SQUAREM
#############################


n <- 10
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 10), 2), round(runif(1, 1, 10)))
  
  epsilon <- 1e-9
  current <- start
  theta <- start
  
  H.current <- -diag(p)
  H.next <- H.current
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    pi <- pi_hat(next1[1], next1[2], next1[3])
    next2 <- update(next1[1], next1[2], next1[3], pi, data)
    
    G.current <- next1 - current
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1 - G.current)
    
    H.next <- (H.current %*% (diag(p) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next
    
    
    theta <- current - H.current %*% G.current
    
    H.current <- dot(u,v)/dot(v,v) * diag(p)
    
    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      iter <- NA
      break
    }
    
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
      count = count+1
    }
    
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    
    current <- theta
  }
  end.time <- Sys.time()
  
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  } else  {
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}

time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)

################################
#### JnJ-QN1
################################


n <- 100
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
fails <- 0
count <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 10), 2), round(runif(1, 1, 10)))
  
  epsilon <- 1e-9
  current <- start
  theta <- start
  
  H.current <- -diag(p)
  H.next <- H.current
  diff <- 100
  iter <- 0
  
  start.time <- Sys.time()
  while(diff > epsilon){
    
    iter <- iter + 1
    print(iter)
    pi <- pi_hat(current[1], current[2], current[3])
    next1 <- update(current[1], current[2], current[3], pi, data)
    
    G.current <- as.matrix(next1 - current)
    theta <- current - H.current %*% G.current
    while(theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      theta <- current + del.theta/2 
      count = count+1
    }
    pi <- pi_hat(theta[1], theta[2], theta[3])
    G.theta <- as.matrix(update(theta[1], theta[2], theta[3], pi, data) - theta)
    
    del.theta <- theta - current
    del.G <- G.theta - G.current
    
    foo <- t(del.theta) %*% H.current
    H.next <- ((del.theta - (H.current %*% del.G)) %*% foo)/as.numeric(foo %*% del.G)
    H.current <- H.next
    
    
    diff <- abs(likelihood(data, theta[1], theta[2], theta[3]) - 
                  likelihood(data, current[1], current[2], current[3]))/
      (abs(likelihood(data, current[1], current[2], current[3])) + 1)
    
    current <- theta
  }
  end.time <- Sys.time()
  if (theta[1]< 0.05 || theta[1] > 0.95 || is.na(iter)){
    fails = fails+1
    time_rep[j] <- NA
    evals[j] <- NA
    estimates[j,] <- NA
  }
  else{
    time_rep [j] <- end.time - start.time
    evals[j] <- iter
    estimates[j,] <- theta
  }
  
}
time_rep <- time_rep[!is.na(time_rep)]
evals <- evals[!is.na(evals)]
estimates <- estimates[complete.cases(estimates), ]
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
mle <- colMeans(estimates)
print(mle)
print(likelihood(data, mle[1], mle[2], mle[3]))
print(fails/n)