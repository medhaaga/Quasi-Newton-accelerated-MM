
log.likelihood <- function(x, batch, freq1, freq2, freq3, freq4){
  
  pi <- x[1]
  alpha <- x[2]
  freq <- freq1 + freq2 + freq3 + freq4
  
  distri <- function(x, batch = 4, pi, alpha){
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
  
  
  denominator <- freq*log(1 - distri(x = 0, batch = 4, pi, alpha))
  dens1 <- freq1*log(distri(x=1, batch=4, pi=pi, alpha))
  dens2 <- freq2*log(distri(x=2, batch=4, pi, alpha))
  dens3 <- freq3*log(distri(x=3, batch=4, pi, alpha))
  dens4 <- freq4*log(distri(x=4, batch=4, pi, alpha))
  
  return (-dens1 - dens2 - dens3 - dens4 + denominator)
}

update <- function(now, batch = 4, freq1, freq2, freq3, freq4){
  
  new <- now
  pi <- now[1]
  alpha <- now[2]
  freq <- freq1 + freq2 + freq3 + freq4
  prod2 <- 1
  prod3 <- 1
  for (k in 0:(batch-1)){
    prod2 <- prod2*(1 - pi + k*alpha)
    prod3 <- prod3*(1 + k*alpha)
  }
  foo <- prod2/prod3
  foo <- foo/(1 - foo)
  s1 <- c(freq, freq2 + freq3 + freq4, freq3 + freq4, freq4)
  s2 <- c(freq1 + freq2 + freq3 + freq*foo, freq1 + freq2 + freq*foo, freq1 + freq*foo, freq*foo)
  r <- rep(freq*(1+foo), 4)
  num1 <- 0
  num2 <- 0
  denom1 <- 0
  denom2 <- 0
  for (k in 0:3){
    num1 <- num1 + (((s1[k+1]*k*alpha)/(pi + k*alpha)) + ((s2[k+1]*k*alpha)/(1 - pi + k*alpha)))
    num2 <- num2 + ((s1[k+1]*pi)/(pi + k*alpha))
    denom1 <- denom1 + ((r[k+1]*k)/(1 + k*alpha))
    denom2 <- denom2 + (((s1[k+1]*pi) / (pi + k*alpha)) + ((s2[k+1]*(1 - pi)) / (1 - pi + k*alpha)))
  }
  new[2] <- num1/denom1
  new[1] <- num2/denom2
  return (new)
}


f <- function(pi, alpha, batch, freq1, freq2, freq3, freq4){
  
  freq <- freq1 + freq2 + freq3 + freq4
  
  distri <- function(x, batch = 4, pi, alpha){
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
  
  
  denominator <- freq*log(1 - distri(x = 0, batch = 4, pi, alpha))
  dens1 <- freq1*log(distri(x=1, batch=4, pi=pi, alpha))
  dens2 <- freq2*log(distri(x=2, batch=4, pi, alpha))
  dens3 <- freq3*log(distri(x=3, batch=4, pi, alpha))
  dens4 <- freq4*log(distri(x=4, batch=4, pi, alpha))
  
  return (-dens1 - dens2 - dens3 - dens4 + denominator)
}

######## function returns TRUE if the parameters satisfy constraint 

param_constraint <- function(par){
  pi <- par[1]
  alpha <- par[2]
  if (alpha > 0 && pi > 0 && pi < 1)
    return(TRUE)
  else
    return(FALSE)
}

