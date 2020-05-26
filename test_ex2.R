set.seed(1)

quadratic <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
    return(result)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
    return(c(x,x))
  }
  else {print("There are no real roots.")} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

rayleigh <- function(x, A, B){
  x <- as.matrix(x)
  num <- t(x) %*% A %*% x
  denom <- t(x) %*% B %*% x
  return(num/denom)
}

update <- function(x, A, B, dir = c("ascent", "descent")){

  u <- as.matrix(x)
  v <- (A - as.numeric(rayleigh(x, A, B))*B) %*% x

  uAu <- t(u) %*% A %*% u
  vAv <- t(v) %*% A %*% v
  uAv <- t(u) %*% A %*% v
  vAu <- t(v) %*% A %*% u
  uBu <- t(u) %*% B %*% u
  vBv <- t(v) %*% B %*% v
  uBv <- t(u) %*% B %*% v
  vBu <- t(v) %*% B %*% u

  a <- (vAv*vBu + vAu*vBv + vAv*uBv) - (vAv*vBu + uAv*vBv + vAu*vBv)
  b <- (vAv*uBu + vAu*vBu + vAu*uBv) - (vAu*vBu + uAv*vBu + uAu*vBv)
  c <- (vAu*uBu) - (uAu*vBu)

  C <- quadratic(a, b, c)
  x1 <- u + C[1]*v
  x2 <- u + C[2]*v
  R1 <- rayleigh(x1, A, B)
  R2 <- rayleigh(x2, A, B)

  if(R1 > R2) {
    ascent <- x1
    descent <- x2
  }
  else {
    ascent <- x2
    descent <- x1
  }

  if (dir == "ascent") {return(ascent)}
  else {return(descent)}

}

update_s <- function(x, A, B, dir = c("ascent", "descent"), s){

  y <- x
  for (i in 1:s){
    y <- update(y, A, B, dir = dir)
  }
  return(y)
}

p <- 10
C <- matrix(runif(p^2, min = -5, max = 5), nrow = p, ncol = p)
D <- matrix(runif(p^2, min = -5, max = 5), nrow = p, ncol = p)
A <- C + t(C)
B <- D %*% t(D)

start <- as.matrix(rep(100, p))
rep <- 100

###########################################
## Naive Algorithm
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
  theta <- update_s(current, A, B, dir = "descent", 2)
  rel.diff <- abs(rayleigh(theta, A, B) - rayleigh(current, A, B))/(abs(rayleigh(current, A, B)) + 1)
  current <- theta
  store_MM <- rbind(store_MM, theta)
}
end.time <- Sys.time()
print(iter)
print(theta)
print(end.time - start.time)

##########################################
## Zhou's quasi-Newton for q=1
##########################################

epsilon <- 1e-9
evals <- rep(0, rep)
time <- rep(0, rep)
store_Zhou <- matrix(0, nrow = rep, ncol = p)

for (i in 1:rep){
  print(i)
  theta <- start
  current <- theta
  H.current <- diag(p)
  H.next <- H.current
  rel.diff <- 100
  iter <- 0

start.time <- Sys.time()
while(rel.diff > epsilon){

  iter <- iter + 1
  next1 <- update_s(current, A, B, dir = "descent", 2)
  next2 <- update_s(next1, A, B, dir = "descent", 2)
  u <- as.matrix(next1 - current)
  v <- as.matrix(next2 - next1)
  foo <- -(dot(u,u))/(dot(u, (v-u)))
  theta <- (1 - foo)*next1 + foo*next2
  if(rayleigh(theta, A, B) < rayleigh(current, A, B)){
    theta <- next2
  }
  rel.diff <- abs(rayleigh(theta, A, B) - rayleigh(current, A, B))/(abs(rayleigh(current, A, B)) + 1)
  current <- theta
}
end.time <- Sys.time()
time[i] <- end.time - start.time
evals[i] <- iter
store_Zhou[i,] <- theta
}

print(mean(time))
print(mean(evals))
print(colMeans(store_BFGS1))




########################################
## Classical BFGS
########################################

epsilon <- 1e-9
evals <- rep(0, rep)
time <- rep(0, rep)
store_BFGS1 <- matrix(0, nrow = rep, ncol = p)

for (i in 1:rep){
  print(i)
  theta <- start
  current <- theta
  H.current <- diag(p)
  H.next <- H.current
  rel.diff <- 100
  iter <- 0

  start.time <- Sys.time()
  while(rel.diff > epsilon){

    iter <- iter + 1
    next1 <- update(current, A, B, dir = "descent")
    next2 <- update(next1, A, B, dir = "descent")

    G.current <- current - next1
    theta <- current - H.current %*% G.current
    if(rayleigh(theta, A, B) < rayleigh(current, A, B)){
      theta <- next2
    }
    G.theta <- theta - update(theta, A, B, dir = "descent")

    u <- as.matrix(theta - current)
    v <- as.matrix(G.theta - G.current)

    H.next <- (H.current %*% (diag(p) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next
    rel.diff <- abs(rayleigh(theta, A, B) - rayleigh(current, A, B))/(abs(rayleigh(current, A, B)) + 1)
    current <- theta
    store_BFGS <- rbind(store_BFGS, as.vector(theta))
}

end.time <- Sys.time()
time[i] <- end.time - start.time
evals[i] <- iter
store_BFGS1[i,] <- theta
}

print(mean(time))
print(mean(evals))
print(colMeans(store_BFGS1))



########################################
## Classical BFGS on Zhou's approach
########################################


epsilon <- 1e-9
evals <- rep(0, rep)
time <- rep(0, rep)
store_BFGS2 <- matrix(0, nrow = rep, ncol = p)
for (i in 1:rep){
  print(i)
  theta <- start
  current <- theta
  H.current <- diag(p)
  H.next <- H.current
  rel.diff <- 100
  iter <- 0

  start.time <- Sys.time()
  while(rel.diff > epsilon){

    iter <- iter + 1
    next1 <- update(current, A, B, dir = "descent")
    next2 <- update(next1, A, B, dir = "descent")
    G.current <- current - next1
    theta <- current - H.current %*% G.current
    if(rayleigh(theta, A, B) < rayleigh(current, A, B)){
      theta <- next2
    }

    u <- as.matrix(next1 - current)
    v <- as.matrix(next1 - next2 - G.current)
    H.next <- (H.current %*% (diag(p) - (v %*% t(v))/dot(v,v))) + (u %*% t(v))/dot(v,v)
    H.current <- H.next
    rel.diff <- abs(rayleigh(theta, A, B) - rayleigh(current, A, B))/(abs(rayleigh(current, A, B)) + 1)
    current <- theta
  }

  end.time <- Sys.time()

  time[i] <- end.time - start.time
  evals[i] <- iter
  store_BFGS2[i,] <- theta
}

print(mean(time))
print(mean(evals))
print(colMeans(store_BFGS2))

