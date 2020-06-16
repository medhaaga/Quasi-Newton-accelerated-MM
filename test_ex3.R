set.seed(1)



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
    pi[i+1,2] <- (p * (mu2 ^ i) * exp(-mu1))/denom[i+1]
  }
  return (pi)
}

update <- function(p, mu1, mu2, pi_hat, data){
  new_p <- dot(data[,2], pi_hat[,1])/sum(data[,2])
  new_mu1 <- dot(data[,1]*data[,2], pi_hat[,1])/dot(data[,2], pi_hat[,1])
  new_mu2 <- dot(data[,1]*data[,2], pi_hat[,2])/dot(data[,2], pi_hat[,2])
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

n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)

for (j in 1:n){

start <- c(round(runif(1, .05, .95), 2), round(runif(1, 1, 100), 2), round(runif(1, 1, 100)))

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
pi <- pi_hat(current[1], current[2], current[3])
theta <- update(current[1], current[2], current[3], pi, data)
diff <- norm((theta-current), type="2")
current <- theta
}
end.time <- Sys.time()
time_rep [j] <- end.time - start.time
evals[j] <- iter
estimates[j,] <- theta
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(colMeans(estimates))


##########################################
## Zhou's quasi-Newton for q=1
##########################################

n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
failures <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .05, .95), 2), round(runif(1, 1, 100), 2), round(runif(1, 1, 100)))

  epsilon <- 1e-7
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
    u <- as.matrix(next1 - current)
    v <- as.matrix(next2 - next1)
    foo <- -(dot(u,u))/(dot(u, (v-u)))
    theta <- (1 - foo)*next1 + foo*next2
    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      failures = failures + 1
      theta <- rep(0, p)
      iter <- 0
      time <- 0
      break
    }

    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
    }
    diff <- norm((theta-current), type="2")
    current <- theta
}
end.time <- Sys.time()
time_rep [j] <- end.time - start.time
evals[j] <- iter
estimates[j,] <- theta
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(colMeans(estimates))




########################################
## Classical BFGS
########################################


n <- 1
time_rep <- rep(0,n)
evals <- rep(0,n)
estimates <- matrix(0, nrow = n, ncol = p)
failures <- 0

for (j in 1:n){
  print(j)
  start <- c(round(runif(1, .05, .95), 2), round(runif(1, 1, 100), 2), round(runif(1, 1, 100)))

  epsilon <- 1e-7
  time <- 0
  count <- 0
  theta <- start
  current <- theta
  H.current <- -diag(p)
  H.next <- H.current
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

    H.next <- (H.current %*% (diag(p) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next


    theta <- current - H.current %*% G.current


    if (theta[1] <= 0 || theta[1] >= 1 || theta[2] <= 0 || theta[3] <= 0){
      print("Failed")
      failures = failures + 1
      theta <- rep(0, p)
      iter <- 0
      time <- 0
      break
    }
    if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
      theta <- next2
      count <- count+1
    }

    diff <- norm((theta-current), type="2")
    current <- theta
  }
  end.time <- Sys.time()
  time_rep [j] <- end.time - start.time
  evals[j] <- iter
  estimates[j,] <- theta
}
print(quantile(time_rep, probs = c(0, .5, 1)))
print(quantile(evals, probs = c(0, .5, 1)))
print(colMeans(estimates))


##########################################
#### SqS1
##########################################



epsilon <- 1e-7
time <- 0
count <- 0
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

  foo <- (dot(u,u))/(dot(u,v))

  theta <- current - 2*foo*u + (foo^2)*v

  if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
    theta <- next2
    count <- count+1
  }

  diff <- norm((theta-current), type="2")
  current <- theta
}

end.time <- Sys.time()
time <- end.time - start.time

print((time))
print((iter))
print(theta)
print(count/iter)


##########################################
#### SqS2
##########################################



epsilon <- 1e-7
time <- 0
count <- 0
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

  if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
    theta <- next2
    count <- count+1
  }

  diff <- norm((theta-current), type="2")
  current <- theta
}

end.time <- Sys.time()
time <- end.time - start.time

print((time))
print((iter))
print(theta)
print(count/iter)

##########################################
#### SqS3
##########################################

epsilon <- 1e-7
time <- 0
count <- 0
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

  foo <- -sqrt((dot(u,u))/(dot(v,v)))

  theta <- current - 2*foo*u + (foo^2)*v

  if(likelihood(data, theta[1], theta[2], theta[3]) < likelihood(data, current[1], current[2], current[3])){
    theta <- next2
    count <- count+1
  }

  diff <- norm((theta-current), type="2")
  current <- theta
}

end.time <- Sys.time()
time <- end.time - start.time

print((time))
print((iter))
print(theta)
print(count/iter)
