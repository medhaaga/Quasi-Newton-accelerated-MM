set.seed(1)
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

result <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

density <- function(x, batch = 4, pi, alpha){
  prod1 <- 1
  prod2 <- 1
  prod3 <- 1
  if (x != 0){
    for (j in 0:(x-1)){
      prod1 <- prod1*(pi + j*alpha)
    }
  }
  if(x != batch){
    for (k in 0:(batch-x-1)){
      prod2 <- prod2*(1 - pi + k*alpha)
    }
  }
  for (l in 0:(batch-1)){
    prod3 <- prod3*(1 + l*alpha)
  }
  return ((choose(batch, x) * prod1 * prod2)/prod3)
}

log.likelihood <- function(freq1, freq2, freq3, freq4, pi, alpha){

  freq <- freq1 + freq2 + freq3 + freq4
  denominator <- freq*log(1 - density(x = 0, batch = 4, pi, alpha))
  dens1 <- freq1*log(density(x=1, batch=4, pi, alpha))
  dens2 <- freq2*log(density(x=2, batch=4, pi, alpha))
  dens3 <- freq3*log(density(x=3, batch=4, pi, alpha))
  dens4 <- freq4*log(density(x=4, batch=4, pi, alpha))

  return (dens1 + dens2 + dens3 + dens4 - denominator)
}

update <- function(batch = 4, current, freq1, freq2, freq3, freq4){
  theta <- current
  foo <- density(x = 0, batch = batch, pi <- current[1], alpha <- current[2])
  foo <- foo/(1 - foo)
  s1 <- c(freq, freq2 + freq3 + freq4, freq3 + freq4, freq4)
  s2 <- c(freq1 + freq2 + freq3 + freq*foo, freq1 + freq2 + freq*foo, freq1 + freq*foo, freq*foo)
  r <- rep(freq*(1+foo), 4)
  num1 <- 0
  num2 <- 0
  denom1 <- 0
  denom2 <- 0
  for (k in 0:3){
    num1 <- num1 + (((s1[k+1]*k*current[2])/(current[1] + k*current[2])) + ((s2[k+1]*k*current[2])/(1 - current[1] + k*current[2])))
    num2 <- num2 + ((s1[k+1]*current[1])/(current[1] + k*current[2]))
    denom1 <- denom1 + ((r[k+1]*k)/(1 + k*current[2]))
    denom2 <- denom2 + (((s1[k+1]*current[1]) / (current[1] + k*current[2])) + ((s2[k+1]*(1 - current[1])) / (1 - current[1] + k*current[2])))
  }
  theta[2] <- num1/denom1
  theta[1] <- num2/denom2
  return (theta)
}

batch <- 4
freq1 <- 26
freq2 <- 15
freq3 <- 3
freq4 <- 9
freq <- freq1 + freq2 + freq3 + freq4
start <- c(.5, 1)

###########################################
## MM Algorithm
###########################################

theta <- start
current <- theta
rel.diff <- 100
epsilon <- 1e-9
iter <- 0
store_MM <- current

start.time <- Sys.time()
while((rel.diff > epsilon))
{
  iter <- iter + 1
  if(iter %% 1000 == 0) print(iter)
  theta <- update(batch, current, freq1, freq2, freq3, freq4)
  rel.diff <- abs(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) - log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))/(abs(log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2])) + 1)
  current <- theta
  store_MM <- rbind(store_MM, theta)
}
end.time <- Sys.time()
print(iter)
print(theta)
print(end.time - start.time)
print(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]))


##########################################
## Zhou's quasi-Newton for q=1
##########################################

theta <- start
current <- theta
rel.diff <- 100
epsilon <- 1e-9
iter <- 0
store_Zhou <- current

start.time <- Sys.time()
while(rel.diff > epsilon){

  iter <- iter + 1
  if (iter%%10 == 0){print(iter)}
  next1 <- update(batch, current, freq1, freq2, freq3, freq4)
  next2 <- update(batch, next1, freq1, freq2, freq3, freq4)
  u <- as.matrix(next1 - current)
  v <- as.matrix(next2 - next1)
  foo <- -(dot(u,u))/(dot(u, (v-u)))
  theta <- (1 - foo)*next1 + foo*next2
  if(theta[1]<=0 || theta[1]>=1 || theta[2]<=0 || (log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) < log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))){
    theta <- next2
  }
  rel.diff <- abs(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) - log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))/(abs(log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2])) + 1)
  current <- theta
  store_Zhou <- rbind(store_Zhou, theta)
}
end.time <- Sys.time()
print(iter)
print(theta)
print(end.time - start.time)
print(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]))


########################################
## Classical BFGS
########################################


theta <- start
current <- theta
H.current <- diag(2)
H.next <- H.current
rel.diff <- 100
epsilon <- 1e-9
iter <- 0
store_BFGS <- current
count = 0

start.time <- Sys.time()
while(rel.diff > epsilon){

  iter <- iter + 1
  if (iter %% 10 == 0){print(iter)}
  next1 <- update(batch, current, freq1, freq2, freq3, freq4)
  G.current <- current - next1
  theta <- current - H.current %*% G.current
  if(round(theta[1], 4) <= 0 || theta[1]>=1 || theta[2]<=0 || (log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) < log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))){
    theta <- update(batch, next1, freq1, freq2, freq3, freq4)
    count = count +1
  }
  G.theta <- theta - update(batch, theta, freq1, freq2, freq3, freq4)

  u <- as.matrix(theta - current)
  v <- as.matrix(G.theta - G.current)
  #H.next <- (diag(2) - (u%*%t(v))/dot(v, u)) %*% H.current %*% (diag(2) - (v %*% t(u))/dot(v, u)) + (u %*% t(u))/dot(v, u)
  H.next <- (H.current %*% (diag(2) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
  H.current <- H.next
  rel.diff <- abs(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) - log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))/(abs(log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2])) + 1)
  current <- theta
  store_BFGS <- rbind(store_BFGS, as.vector(theta))
}

end.time <- Sys.time()
print(theta)
print(iter)
print(H.current)
print(1 - count/iter)
print(end.time - start.time)
print(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]))


########################################
## Classical BFGS like Zhou's
########################################


theta <- start
current <- theta
H.current <- diag(2)
H.next <- H.current
rel.diff <- 100
epsilon <- 1e-9
iter <- 0
store_BFGS <- current
count = 0

start.time <- Sys.time()
while(rel.diff > epsilon){

  iter <- iter + 1
  if (iter %% 10 == 0){print(iter)}
  next1 <- update(batch, current, freq1, freq2, freq3, freq4)
  next2 <- update(batch, next1, freq1, freq2, freq3, freq4)
  G.current <- current - next1
  theta <- current - H.current %*% G.current
  if(round(theta[1], 4) <= 0 || theta[1]>=1 || theta[2]<=0 || (log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) < log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))){
    theta <- next2
    count = count +1
  }
  G.theta <- theta - update(batch, theta, freq1, freq2, freq3, freq4)

  u <- as.matrix(next1 - current)
  v <- as.matrix(next1 - next2 - G.current)

  H.next <- (H.current %*% (diag(2) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
  H.current <- H.next
  rel.diff <- abs(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]) - log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2]))/(abs(log.likelihood(freq1, freq2, freq3, freq4, current[1], current[2])) + 1)
  current <- theta
  store_BFGS <- rbind(store_BFGS, as.vector(theta))
}

end.time <- Sys.time()
print(theta)
print(iter)
print(H.current)
print(1 - count/iter)
print(end.time - start.time)
print(log.likelihood(freq1, freq2, freq3, freq4, theta[1], theta[2]))

