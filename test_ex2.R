library(pracma)
library(SQUAREM)
library(BfgsQN)


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

rayleigh <- function(x, A, B, dir){
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
  num <- t(x1) %*% A %*% x1
  denom <- t(x1) %*% B %*% x1
  R1 <- num/denom
  num <- t(x2) %*% A %*% x2
  denom <- t(x2) %*% B %*% x2
  R2 <- num/denom

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

dim <- 100
C <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
D <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
A <- C + t(C)
B <- D %*% t(D)

start <- as.matrix(rep(1, dim))


###########################################
## Naive Algorithm
###########################################

theta <- start
current <- theta
rel.diff <- 100
epsilon <- 1e-7
feval <- 0
objeval <- 0
dir <- "ascent"

start.time <- Sys.time()
while((rel.diff > epsilon))
{
  feval <- feval + 1
  if(feval %% 1000 == 0) print(feval)
  theta <- update(current, A, B, dir)
  obj.theta <- rayleigh(theta, A, B, dir)
  obj.current <- rayleigh(current, A, B, dir)
  objeval <- objeval + 2
  rel.diff <- abs(obj.theta - obj.current)/(abs(obj.current) + 1)
  current <- theta
}
end.time <- Sys.time()

print(feval)
print(objeval)
print(rayleigh(theta, A, B, dir))
print(end.time - start.time)

##########################################
## Zhou's quasi-Newton for q=1
##########################################

theta <- start
current <- theta
rel.diff <- 100
epsilon <- 1e-5
iter <- 0
feval <- 0
objeval <- 0
fallback <- 0
dir <- "ascent"
max.iter <- 1e4

start.time <- Sys.time()
while(rel.diff > epsilon && iter <= max.iter){

  iter <- iter + 1
  if(iter %% 1000 == 0) print(iter)
  next1 <- update(current, A, B, dir)
  next2 <- update(next1, A, B, dir)
  feval <- feval + 2
  u <- as.matrix(next1 - current)
  v <- as.matrix(next2 - next1)
  alpha <- -(dot(u,u))/(dot(u, (v-u)))
  theta <- (1 - alpha)*next1 + alpha*next2
  if(rayleigh(theta, A, B, dir) < rayleigh(current, A, B, dir)){
    theta <- next2
    fallback <- fallback + 1
    print("fall")
  }
  obj.theta <- rayleigh(theta, A, B, dir)
  obj.current <- rayleigh(current, A, B, dir)
  objeval <- objeval + 2
  rel.diff <- abs(obj.theta - obj.current)/(abs(obj.current) + 1)
  current <- theta
}

end.time <- Sys.time()
print(feval)
print(objeval)
print(rayleigh(theta, A, B, dir))
print(end.time - start.time)
print(fallback/iter)

########################################
### SQUAREM -1
########################################


epsilon <- 1e-7

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K = 1, method=1, maxiter = 1e4, tol = epsilon))
end.time <- Sys.time()

print(fp$value.objfn)
print(fp$fpevals)
print(fp$objfevals)
print(end.time - start.time)

########################################
### SQUAREM -2
########################################


epsilon <- 1e-7

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K = 1, method=2, maxiter = 1e4, tol = epsilon))
end.time <- Sys.time()

print(fp$value.objfn)
print(fp$fpevals)
print(fp$objfevals)
print(end.time - start.time)

########################################
### SQUAREM -3
########################################


epsilon <- 1e-7

start.time <- Sys.time()
fp <- squarem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(K = 1, method=3, maxiter = 1e4, tol = epsilon))
end.time <- Sys.time()

print(fp$value.objfn)
print(fp$fpevals)
print(fp$objfevals)
print(end.time - start.time)

########################################
## Classical BFGS
########################################

epsilon <- 1e-7

start.time <- Sys.time()
fp <- BFGS(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(maxiter = 1e5, tol = epsilon))
end.time <- Sys.time()

print(fp$value.objfn)
print(fp$fpevals)
print(fp$objfevals)
print(end.time - start.time)

##########################################
### L-BFGS
##########################################

epsilon <- 1e-7
dir <- "ascent"

start.time <- Sys.time()
fp <- LBFGS(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir, control = list(m = 10, maxiter = 1e4, tol = epsilon, objfn.inc = .1))
end.time <- Sys.time()

print(fp$value.objfn)
print(fp$fpevals)
print(fp$objfevals)
print(end.time - start.time)

