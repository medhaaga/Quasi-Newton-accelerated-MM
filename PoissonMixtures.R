set.seed(100)
library(pracma)
library(BfgsQN)
library(SQUAREM)


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


likelihood <- function(par, data)
{
  p <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]
  like <- 0
  for (i in 0:9){
    like <- like + data[i+1,2] * log(((p * exp(-mu1) * (mu1 ^ i))/factorial(i)) + (((1-p) * exp(-mu2) * (mu2 ^ i))/factorial(i)))
  }
  return(like)
}


update <- function(par, data){

  p <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]

  pi <- matrix(0, nrow = 10, ncol = 2)
  denom <- rep(0, 10)
  for (i in 0:9){
    denom[i+1] <- (p * (mu1 ^ i) * exp(-mu1)) + ((1-p) * (mu2 ^ i) * exp(-mu2))
    pi[i+1,1] <- (p * (mu1 ^ i) * exp(-mu1))/denom[i+1]
    pi[i+1,2] <- 1 - pi[i+1,1]
  }

  new_p <- dot(data[,2], pi[,1])/sum(data[,2])
  new_mu1 <- dot((data[,1]*data[,2]), pi[,1])/dot(data[,2], pi[,1])
  new_mu2 <- dot((data[,1]*data[,2]), pi[,2])/dot(data[,2], pi[,2])
  return(c(new_p, new_mu1, new_mu2))
}

############################################################3

deaths <- (0:9)
freq <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
data <- as.matrix(cbind(deaths, freq))
truth <- c(.3599, 1.256, 2.663)
#start <- c(.3, 1, 2.5)
dim <- 3

###########################################
## Naive Algorithm
###########################################

N <- 10
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){

start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
start.time <- Sys.time()
fp <- fpiter(par = start, data=data, fixptfn = update, objfn = likelihood, control = list(tol = epsilon, maxiter = 4e3))
end.time <- Sys.time()

if(fp$convergence){
  time_rep [j] <- end.time - start.time
  fpevals[j] <- fp$fpevals
  objevals[j] <- fp$objfevals
  obj.values[j] <- fp$value.objfn
} else{
  time_rep [j] <- NA
  fpevals[j] <- NA
  objevals[j] <- NA
  obj.values[j] <- NA
  fails = fails+1
}
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fpevals, probs = c(0, .5, 1)))
print(quantile(objevals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(fails/N)


##########################################
## Zhou's quasi-Newton for q=1
##########################################

N <- 10
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){

  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))

  epsilon <- 1e-7
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


N <- 10
epsilon <- 1e-9
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
  start.time <- Sys.time()
  fp <- BFGS(par = start, data=data, fixptfn = update, objfn = likelihood, control = list(tol = epsilon, maxiter = 4e3, objfn.inc = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fpevals[j] <- fp$fpevals
    objevals[j] <- fp$objfevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    fpevals[j] <- NA
    objevals[j] <- NA
    obj.values[j] <- NA
    fails = fails+1
  }

}
time_rep <- time_rep[!is.na(time_rep)]
fpevals <- fpevals[!is.na(fpevals)]
objevals <- objevals[!is.na(objevals)]
obj.values <- obj.values[!is.na(obj.values)]

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fpevals, probs = c(0, .5, 1)))
print(quantile(objevals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(fails/N)

##########################################
#### SqS1
##########################################


N <- 10
epsilon <- 1e-9
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
  start.time <- Sys.time()
  fp <- squarem(par = start, data=data, fixptfn = update, objfn = likelihood, control = list(K=1, method = 1, tol = epsilon, maxiter = 4e3, objfn.inc = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fpevals[j] <- fp$fpevals
    objevals[j] <- fp$objfevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    fpevals[j] <- NA
    objevals[j] <- NA
    obj.values[j] <- NA
    fails = fails+1
  }

}
time_rep <- time_rep[!is.na(time_rep)]
fpevals <- fpevals[!is.na(fpevals)]
objevals <- objevals[!is.na(objevals)]
obj.values <- obj.values[!is.na(obj.values)]

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fpevals, probs = c(0, .5, 1)))
print(quantile(objevals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(fails/N)


##########################################
#### SqS2
##########################################



N <- 10
epsilon <- 1e-9
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
  start.time <- Sys.time()
  fp <- squarem(par = start, data=data, fixptfn = update, objfn = likelihood, control = list(K=1, method = 2, tol = epsilon, maxiter = 4e3, objfn.inc = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fpevals[j] <- fp$fpevals
    objevals[j] <- fp$objfevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    fpevals[j] <- NA
    objevals[j] <- NA
    obj.values[j] <- NA
    fails = fails+1
  }

}
time_rep <- time_rep[!is.na(time_rep)]
fpevals <- fpevals[!is.na(fpevals)]
objevals <- objevals[!is.na(objevals)]
obj.values <- obj.values[!is.na(obj.values)]

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fpevals, probs = c(0, .5, 1)))
print(quantile(objevals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(fails/N)


##########################################
#### SqS3
##########################################


N <- 10
epsilon <- 1e-9
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:N){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 50), 2), round(runif(1, 1, 50)))
  start.time <- Sys.time()
  fp <- squarem(par = start, data=data, fixptfn = update, objfn = likelihood, control = list(K=1, method = 3, tol = epsilon, maxiter = 4e3, objfn.inc = 1))
  end.time <- Sys.time()

  if(fp$convergence){
    time_rep [j] <- end.time - start.time
    fpevals[j] <- fp$fpevals
    objevals[j] <- fp$objfevals
    obj.values[j] <- fp$value.objfn
  } else{
    time_rep [j] <- NA
    fpevals[j] <- NA
    objevals[j] <- NA
    obj.values[j] <- NA
    fails = fails+1
  }

}
time_rep <- time_rep[!is.na(time_rep)]
fpevals <- fpevals[!is.na(fpevals)]
objevals <- objevals[!is.na(objevals)]
obj.values <- obj.values[!is.na(obj.values)]

print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(fpevals, probs = c(0, .5, 1)))
print(quantile(objevals, probs = c(0, .5, 1)))
print(quantile(obj.values, probs = c(0, .5, 1)))
print(fails/N)


################################
#### JnJ-QN1
################################


N <- 10
epsilon <- 1e-9
time_rep <- rep(0,N)
fpevals <- rep(0,N)
objevals <- rep(0,N)
obj.values <- rep(0,N)
fails <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .1, .9), 2), round(runif(1, 1, 10), 2), round(runif(1, 1, 10)))

  current <- start
  G.current <- update(current, data) - current
  theta <- start

  H.current <- -diag(dim)
  H.next <- H.current
  diff <- 100
  iter <- 0
  fevals <- 1

  start.time <- Sys.time()
  while(diff > epsilon){

    iter <- iter + 1
    print(iter)
    theta <- current - H.current %*% G.current
    del.theta <- theta - current

    while(theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      theta <- current + del.theta/2
      del.theta <- del.theta/2
      print("ping")
    }

    G.theta <- as.matrix(update(theta, data) - theta)
    fevals <- fevals + 1
    del.G <- G.theta - G.current

    foo <- t(del.theta) %*% H.current
    H.next <- H.current + ((del.theta - (H.current %*% del.G)) %*% foo)/as.numeric(foo %*% del.G)
    H.current <- H.next

    diff <- norm(G.theta, type = "2")
    current <- theta
    G.current <- G.theta
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
