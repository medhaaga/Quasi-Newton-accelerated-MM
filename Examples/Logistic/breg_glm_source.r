dg <- cmpfun(function(x, xa, v, plist1, y, h, hgrad, dd) {
  n <- length(plist1)
  p <- length(x)
  u <- double(p)
  
  c <- as.numeric( dd(h(x)) ) 
  # ddGhx <- diag(as.numeric(h(x)))
  for (i in 1:n) {
    z <- as.matrix(plist1[[i]](xa))
    u <- u + v[i]*z
  }
  sum1 <- sum(v)*x - u
  sum2 <- c*( h(x)-y )
  return(sum1 + hgrad(x)%*%sum2) # 
})

ddg <- cmpfun(function(x, v, h, hgrad, dd) {
  p <- length(x)
  u <- matrix(0, p, p)
  dh <- as.matrix(hgrad(x))
  c <- as.numeric( dd( h(x) ) )
  return( sum(v)*diag(p) + (dh*c) %*% t(dh))
})

#Note: make D vectorized
proximity = cmpfun(function(x,v,plist1,y,breg,h) {
  n <- length(plist1)
  p <- length(x)
  obj <- 0
  for (i in 1:n) {
    z <- as.matrix(x - plist1[[i]](x))
    obj <- obj + v[i]*norm(z,'f')**2
  }
  #obj <- obj + mean( D(y,h(x)) )
  obj <- obj + sum( breg(y,x ) )
  return(obj)
}) 

#NEED TO FIX WOOD INV
#following are from package splitfeas
wood_inv_solve = cmpfun(function(x,v,hgrad,df,w=1) {
  n <- length(x)
  vv <- sum(v); ww <- sum(w)
  dh <- as.matrix(hgrad(x))
  p <- dim(dh)[2]  
  return( (1/vv)*df - (ww/vv^2)*dh %*% solve(diag(p) + (ww/vv)*t(dh) %*% dh , t(dh)%*%df ) )
})


backtrack = cmpfun(function(x,dx,f,df,alpha=0.05,beta=0.5) {
  t <- 1
  g <- df(x)
  u <- alpha*sum(g*dx)
  k <- 1
  repeat {
    if (f(x + t*dx) <= f(x) + t*u) break
    t <- beta*t
    print(paste0("backtrack ",k))
    k <- k + 1
  }
  return(t)
})

qnamm = cmpfun(function(x, fx_mm, qn, fx_obj, max_iter=50, tol=1e-6, ...) {
  n = length(x)
  U = matrix(0,n,qn)
  V = matrix(0,n,qn)
  objval = Inf
  objective = double(max_iter)
  objective[1] = fx_obj(x, ...)
  #  Xhist = vector(mode="list", length=qn+max_iter)
  Xhist <- matrix(NA,n,qn+max_iter)
  #
  #   accumulate the first QN differences for Quasi-Newton acceleration  
  #  
  for (i in 1:qn) {
    Xhist[,i] = x
    #    Xhist[[i]] = x
    x_old = x    
    x = fx_mm(x, ...)
    U[,i] = x - x_old
  }
  V[,1:(qn-1)] = U[,2:qn]
  x_old = x    
  x = fx_mm(x, ...)
  V[,qn] = x - x_old
  old_secant = 1
  C = t(U)%*%(U-V)
  nacc = 0
  nrej = 0
  for (i in 1:max_iter) {
    Xhist[,qn+i] = x
    objval_old = objval
    x_old = x
    x = fx_mm(x, ...)      
    #
    #   do one more MM step to accumulate secant pairs  
    #
    U[,old_secant] = x - x_old
    x_old = x
    x = fx_mm(x, ...)
    V[,old_secant] = x - x_old
    C[old_secant,] = t(U[,old_secant,drop=FALSE]) %*% (U-V)
    C[,old_secant] = t(U) %*% (U[,old_secant,drop=FALSE] - V[,old_secant,drop=FALSE])
    new_secant = old_secant
    old_secant = (old_secant %% qn) + 1
    objval_MM = fx_obj(x, ...)   
    #  
    #   quasi-Newton jump
    #     
    #      x_qn = x_old + V %*% solve(C, t(U)%*%U[,new_secant,drop=FALSE])
    x_qn = x_old + V %*% pseudoinverse(C) %*% (t(U)%*%U[,new_secant,drop=FALSE])
    x_qn = fx_mm(x_qn, ...)
    objval_QN = fx_obj(x_qn, ...)
    #
    #     choose MM vs QN jump
    #    
    if (objval_QN < objval_MM) {
      x = x_qn;
      objval = objval_QN;
      nacc = nacc + 1
      #        print('QN step accepted.')
    } else {
      objval = objval_MM;
      #        print('QN step rejected.')
      nrej = nrej + 1
    }
    objective[i+1] = objval
    #    
    # stopping rule
    #    
    print(norm(as.matrix(x-x_old),'f'))
    if (norm(as.matrix(x-x_old),'f') < tol) break
  }
  print(paste("Accepted:", nacc))
  print(paste("Rejected:", nrej))
  return(list(x=x, objective=objective[1:i+1], iter=i+qn, Xhist=Xhist[,1:(i+qn),drop=FALSE]))
})

mmqn_step = cmpfun(function(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=FALSE) {
  df <- dg(x,x,v,plist1,y,h,hgrad,dd)
  if(woodbury==TRUE){
    #    H <- wood_inv(x,v,w,hgrad)
    #    dx <- -H %*% df
    dx <- -wood_inv_solve(x,v,hgrad,df)
  }else{
    H <- ddg(x,v,h,hgrad,dd)
    dx <- -solve(H,df)
  }
  t <- backtrack(x,dx,f,df)
  return(x + t*dx)
})

######### MM Algorithm #############

nmsfp_mm <- cmpfun(function(x0,v,plist1,y,f,df,h,hgrad,dd,woodbury=TRUE,tol=1e-8,max_iter=1e3) {
  x_old <- x0
  x_new <- x0
  diff <- 10
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  iter <- 1
  while(diff > tol && iter < max_iter){
    print(iter)
    x_new <- fx_mm(x_old)
    diff <- norm(x_new - x_old, '2')
    x_old <- x_new
    print(paste('Objective at iter', iter, ' is ', f(x_new), 'diff is ', diff))
    iter <- iter + 1
  }
  return(list(x = x_new, loss = f(x_new), fpevals = iter))
})

########### MM-QN acceleration ###########
nmsfp_mmqn <- cmpfun(function(x0,v,plist1,y,f,df,h,hgrad,dd,woodbury=TRUE,qn=5,tol=1e-8,max_iter=1e3) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  sol <- qnamm(x, fx_mm, qn, f, max_iter=max_iter, tol=tol) 
  return(list(x=sol$x, loss=sol$objective, fpevals=sol$iter))
})

########### BQN acceleration ###########
nmsfp_bqn <- cmpfun(function(x0, v, plist1, y, f, df, h, hgrad, dd, woodbury=TRUE, qn=5, tol=1e-8, max_iter=1e3, step.min=1e-5, objfn.inc = .01, intermed=TRUE) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  sol <- BQN(par = x, fixptfn = fx_mm, objfn = f, 
             control = list(qn=qn, step.min=step.min, tol = tol, maxiter = max_iter, objfn.inc = objfn.inc, verbose=TRUE, intermed=intermed))

  return(list(x=sol$par,loss=sol$value.objfn, fevals=sol$fpevals, Xhist=sol$p.inter))
})

########### LBQN acceleration ###########
nmsfp_lbqn <- cmpfun(function(x0, v, plist1, y, f, df, h, hgrad, dd, woodbury=TRUE, m=10, tol=1e-8, max_iter=1e3, step.min=1e-5, objfn.inc = .01, intermed=TRUE) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  sol <- LBQN(par = x, fixptfn = fx_mm, objfn = f, 
             control = list(m=m, tol = tol, maxiter = max_iter, objfn.inc = objfn.inc, verbose=TRUE, intermed=intermed))
  
  return(list(x=sol$par,loss=sol$value.objfn, fevals=sol$fpevals, Xhist=sol$p.inter))
})

########### SQUAREM acceleration ###########
nmsfp_sq <- cmpfun(function(x0, v, plist1, y, f, df, h, hgrad, dd, woodbury=TRUE, method=1, tol=1e-8, max_iter=1e3, intermed=TRUE) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  sol <- squarem(par = x, fixptfn = fx_mm, objfn = f,
             control = list(K=1, method = method, tol = tol, maxiter = max_iter, trace=TRUE, intermed=intermed))
  return(list(x=sol$par, loss=sol$value.objfn, fevals=sol$fpevals, Xhist=sol$p.intermed))
})

######### DAAREM acceleration ##########
nmsfp_dar <- cmpfun(function(x0, v, plist1, y, f, df, h, hgrad, dd, woodbury=TRUE, tol=1e-8, max_iter=1e3) {
  x <- x0
  fx_mm <- function(x) {return(mmqn_step(x,v,plist1,y,f,df,h,hgrad,dd,woodbury=woodbury))}
  dar_obj <- function(x) {return(-f(x))}
  sol <- daarem(par = x, fixptfn = fx_mm, objfn = dar_obj, control = list(tol = tol, maxiter = max_iter, intermed=TRUE))
  return(list(x=sol$par, loss=sol$value.objfn, fevals=sol$fpeval, Xhist=sol$p.intermed))
})


#project onto the k-sparsity set of x: takes k largest components and rest are zeros
project_sparsity <- function(x,spar){
  res <- rep(0, length(x))
  ind_largest <-  order(abs(x),decreasing=T)[1:spar] #check whether absolute value should be here
  res[ ind_largest ] <- x[ind_largest]
  return(res)
}
