#################################################
######## Generalised Eigenvalues ################
#################################################

rm(list = ls())

rayleigh <- function(x, A, B, dir){
  x <- as.matrix(x)
  num <- t(x) %*% A %*% x
  denom <- t(x) %*% B %*% x
  if (dir == "descent")
    return(num/denom) else
      return(-num/denom)
}

neg.objective <- function(x, A, B, dir){
  x <- as.matrix(x)
  num <- t(x) %*% A %*% x
  denom <- t(x) %*% B %*% x
  if (dir == "descent")
    return(-num/denom) else
      return(num/denom)
}

daarem.objective <- function(x, A, B, dir){
  x <- as.matrix(x)
  num <- t(x) %*% A %*% x
  denom <- t(x) %*% B %*% x
  if (dir == "descent")
    return(num/denom) else
      return(-num/denom)
}

update <- function(x, A, B, dir = c("ascent", "descent")){
  x <- as.matrix(x)
  u <- as.matrix(x)
  v <- (A - as.numeric((t(x) %*% A %*% x)/(t(x) %*% B %*% x))*B) %*% x

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

  delta <- (b^2 - 4*a*c)
  if(delta > 0){ # first case D>0
    x_1 = (-b+sqrt(delta))/(2*a)
    x_2 = (-b-sqrt(delta))/(2*a)
    C <- c(x_1,x_2)
  }
  else if(delta == 0){ # second case D=0
    x = -b/(2*a)
    C <- (c(x,x))
  }
  else {stop("No roots")} # third case D<0

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

