#' Quasi-Newton acceleration of MM algorithm
#'
#' \code{qnamm} performs Quasi-Newton acceleration of an MM algorithm.
#'
#' @param x initial iterate
#' @param fx_mm MM algorithm map
#' @param qn number of secants
#' @param fx_obj handle to objective function
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param V Additional arguments to pass to \code{fx_mm}
#' @import compiler corpcor
#' @export
#' @references H Zhou, D Alexander, and K Lange. (2011) A quasi-Newton acceleration method for high-dimensional optimization algorithms, Statistics and Computing, 21(2):261-273.

library("corpcor")
qnamm <- function(x, fx_mm, qn, fx_obj, max_iter=100, tol=1e-6, V=V) {

  conv <- TRUE

  n <- length(x)
  U <- matrix(0,n,qn)
  W <- matrix(0,n,qn)
  objval <- fx_obj(x, V=V)
  objective <- double(max_iter)
  Xhist <- matrix(NA,n+1,qn+max_iter)
  #
  #   accumulate the first QN differences for Quasi-Newton acceleration
  #
  for (i in 1:qn) {
    Xhist[,i] <- c(x, objval)
    x_old <- x
    x <- fx_mm(x, V=V)
    objval <- fx_obj(x, V=V)
    U[,i] <- x - x_old
  }

  if(qn>1)
    W[,1:(qn-1)] <- U[,2:qn]
  x_old <- x
  x <- fx_mm(x, V=V)
  objval <- fx_obj(x, V=V)
  W[,qn] <- x - x_old

  fevals <- qn+1


  old_secant <- 1
  C <- t(U)%*%(U-W)
  nacc <- 0
  nrej <- 0

  for (i in 1:max_iter) {

    Xhist[,qn+i] <- c(x, objval)
    objval_old <- objval
    x_old <- x
    x <- fx_mm(x, V=V)

    #
    #   do one more MM step to accumulate secant pairs
    #

    U[,old_secant] <- x - x_old
    x_old <- x
    x <- fx_mm(x, V=V)
    W[,old_secant] <- x - x_old
    C[old_secant,] <- t(U[,old_secant,drop=FALSE]) %*% (U-W)
    C[,old_secant] <- t(U) %*% (U[,old_secant,drop=FALSE] - W[,old_secant,drop=FALSE])
    new_secant <- old_secant
    old_secant <- (old_secant %% qn) + 1
    objval_MM <- fx_obj(x, V=V)
    #
    #   quasi-Newton jump
    #
    #      x_qn <- x_old + V %*% solve(C, t(U)%*%U[,new_secant,drop=FALSE])
    x_qn <- x_old + W %*% pseudoinverse(C) %*% (t(U)%*%U[,new_secant,drop=FALSE])
    x_qn <- fx_mm(x_qn, V=V)
    objval_QN <- fx_obj(x_qn, V=V)
    #
    #     choose MM vs QN jump
    #
    if (objval_QN < objval_MM) {
      x <- x_qn;
      objval <- objval_QN;
      nacc <- nacc + 1
    } else {
      objval <- objval_MM;
      nrej <- nrej + 1
    }
    objective[i] <- objval
    #
    # stopping rule
    #
    #print(norm(as.matrix(x-x_old),'f')/(norm(as.matrix(x_old),'f')+1))

    if (norm(as.matrix(x-x_old),"2") < tol) break
  }
  fevals <- fevals + 3*i
  levals <- 2*i
  if (i == max_iter) conv = FALSE
  print(paste("Accepted:", nacc))
  print(paste("Rejected:", nrej))
  return(list(convergence = conv, x=x, fevals = fevals, levals = levals,  accept = nacc, reject = nrej, objective=objective[i], Xhist=Xhist[,1:(i+qn),drop=FALSE]))
}
