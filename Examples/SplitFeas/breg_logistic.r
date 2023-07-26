rm(list=ls())
library(compiler)
library(corpcor)
library(turboEM)
library(daarem)
library(SQUAREM)
library(quasiNewtonMM)
source('breg_glm_source.R')
set.seed(123)

#Logistic Breg GLM
expit = function(x){
  return( exp(x)/(1+exp(x)) )
}

plist1 <- vector(mode = "list", 1)
plist1[[1]] <- function(x, spar) {return(project_sparsity(x, spar=k))}

#CHeck H and DH!!!
#here h is a linear function
h <- function(x) {
  return( expit(A %*% x) )
}

#check this gradient
dh <- function(x){ 
  return( t( A) %*% diag( as.vector( exp(A%*%x)/( exp(A%*%x)+1 )^2 )) )  
  #return( t( diag(as.vector(exp(A%*%x))) %*% A  ) )
  
  #return( t(A)*exp( t(A) * as.vector(x) )  )
}

D2 <- function(y,x){
  aa <- A%*%x
  return( log(1+exp(aa)) - y*aa) #  + y*log(y) - y )
}

#this should provide the diagonal of the hessian, assuming that the hessian has diagonal form!
dd <- function(z,N=1){
  return(N/(N*z-z^2) )
}

# goal function
f <- function(x) {
  return(proximity(x, v, plist1,y, D2, h))
}
# df
df <- function(x) {return(dg(x, x, v, plist1, y, h, dh, dd))}


###### Initialize data and solutions

#overdetermined example
m <- 1500  #rows
n <- 2000 #columns
k = 10    #sparsity
N <- 10

#create true solution
x <- rep(0,n)
x[sample(n,k)] <- (runif(k,4,5))    # 5#*sign(runif(k,-1,1))

#design matrix
A <- matrix(rnorm(n*m,0,.2),m,n) 
pr = expit(A%*%x)
y = rbinom(m,1,pr)
v <- 1e2 #this is the weight in front of the sparsity projection term
x0 <- 1.5*matrix(runif(N*n), N, n)
max_iter <- 1e3
tol = 1e-6
nStarts <- nrow(x0)
woodbury = TRUE
################ Benchmarks ################

#MCP
library(ncvreg)
mcpfit <- cv.ncvreg(A, y, family='binomial')
x_mcp <- coef(mcpfit)[-1]
#SCAD
scadfit <- cv.ncvreg(A, y, family='binomial', penalty = 'SCAD')
x_scad <- coef(scadfit)[-1]
#LASSO
library(glmnet)
cvfit <- glmnet::cv.glmnet(A, y, family = 'poisson')
x_lasso <-coef(cvfit,s='lambda.min',exact=TRUE)[-1] #drop the intercept term

#############################################
################ MM-QN  ####################
#############################################


time_mm <- rep(NA, nStarts)
obj_mm <- rep(NA, nStarts)
eval_mm <- rep(NA,nStarts)
sol_mm <- matrix(NA, n, nStarts)
objlist_mm <- NA

for (i in 1:nStarts){
  time_mm[i] <- system.time(fp <- nmsfp_mm(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, tol = tol))[1]
  sol_mm[,i] <- fp$x
  obj_mm[i] <- fp$loss
  eval_mm[i] <- fp$fpevals
}
time_mm <- time_mm/60

print(round(quantile(eval_mm[!is.na(eval_mm)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_mm[!is.na(time_mm)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_mm[!is.na(obj_mm)], c(.5, 0.25, 0.75)), 3))

save(time_mm, obj_mm, eval_mm, sol_mm, objlist_mm, file = "Out/objects_mm.Rdata")

#############################################
#################### ZAL  ##################
#############################################

sol_feas <- matrix(NA, n, nStarts)
obj_feas <- rep(NA, nStarts)
time_feas <- rep(NA, nStarts)
eval_feas <- rep(NA, nStarts)

for ( i in 1:nStarts) {
  print(i)
  start <- Sys.time()
  solution <- try(nmsfp_mmqn(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, tol = tol), silent=TRUE)
  if (!inherits(solution, "try-error"))
  {
    time_feas[i] <- Sys.time() - start
    sol_feas[, i] <- solution$x[, ncol(solution$x)]  
    obj_feas[i] <- solution$loss[length(solution$loss)]
    eval_feas[i] <- solution$fpevals
    objlist_feas <- solution$loss
  }
}
save(time_feas, sol_feas, obj_feas, objlist_feas, eval_feas, file = "Out/objects_feas.Rdata")
print(round(quantile(eval_feas[!is.na(eval_feas)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_feas[!is.na(time_feas)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_feas[!is.na(obj_feas)], c(.5, 0.25, 0.75)), 3))

#############################################
############ BQN, q = 2 ##############
#############################################


time_bqn1 <- rep(NA, nStarts)
obj_bqn1 <- rep(NA, nStarts)
eval_bqn1 <- rep(NA,nStarts)
sol_bqn1 <- matrix(NA, n, nStarts)

for (i in 1:nStarts){
  print(i)
  start <- Sys.time()
  fp <- nmsfp_bqn(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury=woodbury, qn=2, 
                  tol=tol, max_iter=max_iter, step.min=1, objfn.inc = .01, intermed=TRUE)
  time_bqn1[i] <- Sys.time() - start
  sol_bqn1[,i] <- fp$x
  obj_bqn1[i] <- fp$loss
  eval_bqn1[i] <- fp$fevals
  objlist_bqn1 <- fp$Xhist
}
save(time_bqn1, sol_bqn1, obj_bqn1, objlist_bqn1, eval_bqn1, file = "Out/objects_bqn1.Rdata")
print(quantile(round(eval_bqn1[!is.na(eval_bqn1)],1), c(.25,.5, 0.75)))
print(quantile(round(time_bqn1[!is.na(time_bqn1)],2), c(.25,.5, 0.75)))
print(quantile(round(obj_bqn1[!is.na(obj_bqn1)], 1), c(.25,.5, 0.75)))

#############################################
############ LBQN, m=5  ##############
#############################################


time_lbqn <- rep(NA, nStarts)
obj_lbqn <- rep(NA, nStarts)
eval_lbqn <- rep(NA,nStarts)
sol_lbqn <- matrix(NA, n, nStarts)

for (i in 1:nStarts){
  print(i)
  start <- Sys.time()
  fp <- try(nmsfp_lbqn(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury=woodbury, 
                   m=5, tol=tol, max_iter=max_iter, objfn.inc = .01, intermed=TRUE), silent=TRUE)
  if (!inherits(fp, "try-error"))
  {
    time_lbqn[i] <- Sys.time() - start
    sol_lbqn[,i] <- fp$x
    obj_lbqn[i] <- fp$loss
    objlist_lbqn <- fp$Xhist
    eval_lbqn[i] <- fp$fevals
  }

}
save(time_lbqn, sol_lbqn, obj_lbqn, objlist_lbqn, eval_lbqn, file = "Out/objects_lbqn.Rdata")
print(round(quantile(eval_lbqn[!is.na(eval_lbqn)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_lbqn[!is.na(time_lbqn)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_lbqn[!is.na(obj_lbqn)], c(.5, 0.25, 0.75)), 3))

#############################################
############ LBQN, m=10  ##############
#############################################


time_lbqn2 <- rep(NA, nStarts)
obj_lbqn2 <- rep(NA, nStarts)
eval_lbqn2 <- rep(NA,nStarts)
sol_lbqn2 <- matrix(NA, n, nStarts)

for (i in 1:nStarts){
  start <- Sys.time()
  fp <- nmsfp_lbqn(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury=woodbury, 
                   m=10, tol=tol, max_iter=max_iter, objfn.inc = .01, intermed=TRUE)
  time_lbqn2[i] <- Sys.time() - start
  sol_lbqn2[,i] <- fp$x
  obj_lbqn2[i] <- fp$loss
  objlist_lbqn2 <- fp$Xhist
  eval_lbqn2[i] <- fp$fevals
}
save(time_lbqn2, sol_lbqn2, obj_lbqn2, objlist_lbqn2, eval_lbqn2, file = "Out/objects_lbqn2.Rdata")
print(round(quantile(eval_lbqn2[!is.na(eval_lbqn2)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_lbqn2[!is.na(time_lbqn2)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_lbqn2[!is.na(obj_lbqn2)], c(.5, 0.25, 0.75)), 3))

#############################################
############ LBQN, m=2  ##############
#############################################


time_lbqn1 <- rep(NA, nStarts)
obj_lbqn1 <- rep(NA, nStarts)
eval_lbqn1 <- rep(NA,nStarts)
sol_lbqn1 <- matrix(NA, n, nStarts)

for (i in 1:nStarts){
  start <- Sys.time()
  fp <- nmsfp_lbqn(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury=woodbury, 
                   m=2, tol=tol, max_iter=max_iter, objfn.inc = .01, intermed=TRUE)
  time_lbqn1[i] <- Sys.time() - start
  sol_lbqn1[,i] <- fp$x
  obj_lbqn1[i] <- fp$loss
  objlist_lbqn1 <- fp$Xhist
  eval_lbqn1[i] <- fp$fevals
}
save(time_lbqn1, sol_lbqn1, obj_lbqn1, objlist_lbqn1, eval_lbqn1, file = "Out/objects_lbqn1.Rdata")


#############################################
################ SQUAREM - 1  ###############
#############################################


time_sq1 <- rep(NA, nStarts)
obj_sq1 <- rep(NA, nStarts)
eval_sq1 <- rep(NA,nStarts)
sol_sq1 <- matrix(NA, n, nStarts)
objlist_sq1 <- NA

for (i in 1:nStarts){
  start <- Sys.time()
  fp <- nmsfp_sq(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, method=1, tol = tol, max_iter = max_iter, intermed=TRUE)
  time_sq1[i] <- Sys.time() - start
  obj_sq1[i] <- fp$loss
  eval_sq1[i] <- fp$fevals
  objlist_sq1 <- fp$Xhist
}
save(time_sq1, sol_sq1, obj_sq1, objlist_sq1, eval_sq1, file = "Out/objects_sq1.Rdata")
print(quantile(round(eval_sq1[!is.na(eval_sq1)],0), c(.5, 0.25, 0.75)))
print(quantile(round(time_sq1[!is.na(time_sq1)],2), c(.5, 0.25, 0.75)))
print(quantile(round(obj_sq1[!is.na(obj_sq1)], 2), c(.5, 0.25, 0.75)))


#############################################
################ SQUAREM - 2  ###############
#############################################


time_sq2 <- rep(NA, nStarts)
obj_sq2 <- rep(NA, nStarts)
eval_sq2 <- rep(NA,nStarts)
sol_sq2 <- matrix(NA, n, nStarts)
objlist_sq2 <- NA

for (i in 1:nStarts){
  start <- Sys.time()
  fp <- nmsfp_sq(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, method=2, tol = tol, max_iter = max_iter, intermed=TRUE)
  time_sq2[i] <- Sys.time() - start
  obj_sq2[i] <- fp$loss
  eval_sq2[i] <- fp$fevals
  objlist_sq2 <- fp$Xhist
}
save(time_sq2, sol_sq2, obj_sq2, objlist_sq2, eval_sq2, file = "Out/objects_sq2.Rdata")
print(round(quantile(eval_sq2[!is.na(eval_sq2)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_sq2[!is.na(time_sq2)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_sq2[!is.na(obj_sq2)], c(.5, 0.25, 0.75)), 3))


#############################################
################ SQUAREM - 3  ###############
#############################################


time_sq3 <- rep(NA, nStarts)
obj_sq3 <- rep(NA, nStarts)
eval_sq3 <- rep(NA,nStarts)
sol_sq3 <- matrix(NA, n, nStarts)
objlist_sq3 <- NA

for (i in 1:nStarts){
  start <- Sys.time()
  fp <- nmsfp_sq(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, method=3, tol = tol, max_iter = max_iter, intermed=TRUE)
  time_sq3[i] <- Sys.time() - start
  obj_sq3[i] <- fp$loss
  eval_sq3[i] <- fp$fevals
  objlist_sq3 <- fp$Xhist
}
save(time_sq3, sol_sq3, obj_sq3, objlist_sq3, eval_sq3, file = "Out/objects_sq3.Rdata")
print(round(quantile(eval_sq3[!is.na(eval_sq3)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_sq3[!is.na(time_sq3)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_sq3[!is.na(obj_sq3)], c(.5, 0.25, 0.75)), 3))

#############################################
################ DAAREM  ####################
#############################################


time_dar <- rep(NA, nStarts)
obj_dar <- rep(NA, nStarts)
eval_dar <- rep(NA,nStarts)
sol_dar <- matrix(NA, n, nStarts)

for (i in 1:nStarts){
  print(i)
  start <- Sys.time()
  fp <- try(nmsfp_dar(x0[i,], v, plist1, y, f, df, h, dh, dd, woodbury = woodbury, tol = tol, max_iter = max_iter))
  if (!inherits(fp, "try-error")){
    time_dar[i] <- Sys.time() - start
    obj_dar[i] <- -fp$loss
    objlist_dar <- fp$Xhist
    eval_dar[i] <- fp$fevals
    sol_dar[,i] <- fp$x
  }

}
save(time_dar, sol_dar, obj_dar, objlist_dar, eval_dar, file = "Out/objects_dar.Rdata")
print(round(quantile(eval_dar[!is.na(eval_dar)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(time_dar[!is.na(time_dar)], c(.5, 0.25, 0.75)), 3))
print(round(quantile(obj_dar[!is.na(obj_dar)], c(.5, 0.25, 0.75)), 3))

save(time_mm, time_bqn1, time_lbqn, time_lbqn2, time_sq1, time_sq2, time_sq3, time_dar, time_feas,
     obj_mm, obj_bqn1, obj_lbqn, obj_lbqn2, obj_sq1, obj_sq2, obj_sq3, obj_dar, obj_feas,
     eval_mm, eval_bqn1, eval_lbqn, eval_lbqn2, eval_sq1, eval_sq2, eval_sq3, eval_dar, eval_feas,
     sol_mm, sol_bqn1, sol_lbqn, sol_lbqn2, sol_sq1, sol_sq2, sol_sq3, sol_dar, sol_feas,
     objlist_mm, objlist_bqn1, objlist_lbqn, objlist_lbqn2, objlist_sq1, objlist_sq2, objlist_sq3, objlist_dar, objlist_feas, file = "Out/exp_objects.Rdata")


load(file = "Out/exp_objects.Rdata")

y_range <- range(log(objlist_bqn1[-(1:3),-(1:ncol(objlist_bqn1)-1)]), log(objlist_lbqn[-(1:3),-(1:ncol(objlist_lbqn)-1)]),
                 log(objlist_lbqn2[-(1:3),-(1:ncol(objlist_lbqn2)-1)]), log(objlist_feas[-(1:3)]),
                 log(objlist_sq1[-(1:3),-(1:ncol(objlist_sq1)-1)]), log(-objlist_dar[-(1:3),1]))
pdf(file = "Out/split_feas-objVSiter.pdf", height=5, width=7)
plot(4:nrow(objlist_bqn1), log(objlist_bqn1[-(1:3),-(1:ncol(objlist_bqn1)-1)]), type = "l",
     lwd=3, xlab="Iterations", ylab="Objective Value", col="darkolivegreen3", lty=1, ylim = y_range)
lines(4:(nrow(objlist_lbqn)), log(objlist_lbqn[-(1:3),-(1:ncol(objlist_lbqn)-1)]), col="darkgoldenrod1", lwd=3, lty=2)
lines(4:(nrow(objlist_lbqn2)), log(objlist_lbqn2[-(1:3),-(1:ncol(objlist_lbqn2)-1)]), col="aquamarine3", lwd=3, lty=3)
lines(4:(length(objlist_feas)), log(objlist_feas[-(1:3)]), col="cadetblue1", lwd=3, lty=4)
lines(4:(nrow(objlist_sq1)), log(objlist_sq1[-(1:3),-(1:ncol(objlist_sq1)-1)]), col="darkgreen", lwd=3, lty=5)
lines(4:(nrow(objlist_dar)), log(-objlist_dar[-(1:3),1]), col="brown1", lwd=3, lty=6)
legend("topright", legend=c("BQN, q=2", "LBQN, m=5", "LBQN, m=10", "ZAL", "SqS1", "DAAREM"), 
       col=c("darkolivegreen3", "darkgoldenrod1", "aquamarine3", "cadetblue1", "darkgreen", "brown1"), lwd=3, lty=1:6)
dev.off()


pdf(file = "Out/split_feas-timeVSexp.pdf", height=5, width=7)
yrange = range(time_bqn1, time_lbqn, time_sq1, time_sq2, time_sq3, time_dar, time_feas[!is.na(time_feas)])

plot((1:N), time_bqn1, type='o', lwd=2, col="aquamarine3", ylim=yrange, xlab = "Experiment", ylab="Time (in min)", pch=1)
lines((1:N), time_lbqn, type='o', lwd=2, col="coral", pch=2)
lines((1:N), time_sq1, type='o', lwd=2, col="cadetblue1", pch=3)
lines((1:N), time_dar, type='o', lwd=2, col="darkgoldenrod1", pch=4)
lines((1:N), time_feas, type='o', lwd=2, col="darkgreen", pch=5)
legend("topright", pch=(1:5), legend = c("BQN", "L-BQN", "SQ1", "DAAREM", "ZAL"),
       col = c("aquamarine3", "coral", "cadetblue1", "darkgoldenrod1", "darkgreen"), lwd=2)
dev.off()
pdf(file = "Out/split_feas-objVSexp.pdf", height=5, width=7)
yrange = range(obj_mm, obj_bqn1, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_dar)
plot((1:N), obj_mm, type='o', lwd=2, col="darkolivegreen3", ylim=yrange, pch=6, xlab = "Experiment", ylab = "Objective value")
lines((1:N), obj_bqn1, type='o', lwd=2, col="aquamarine3", pch=1)
lines((1:N), obj_lbqn, type='o', lwd=2, col="coral", pch=2)
lines((1:N), obj_sq1, type='o', lwd=2, col="cadetblue1", pch=3)
lines((1:N), obj_dar, type='o', lwd=2, col="darkgoldenrod1", pch=4)
lines((1:N), obj_feas, type='o', lwd=2, col="darkgreen", pch=5)
lines((1:N), obj_mm, type='o', lwd=2, col="darkolivegreen3", pch=6)
legend("topleft", pch=c(6,(1:5)), legend = c("MM", "BQN", "L-BQN", "SQ1", "DAAREM", "ZAL"),
       col = c("darkolivegreen3", "aquamarine3", "coral", "cadetblue1", "darkgoldenrod1", "darkgreen"), lwd=2)
dev.off()


#signal
round ( sort( abs(x)[abs(x)>0],decreasing=T), 4)
#comparison
sort(round(abs(sol_feas[,1]),4),decreasing=T)[1:(k+2)]
sort(round(abs(sol_bqn1[,1]),4),decreasing=T)[1:(k+2)]
round( sort(abs(x_lasso),decreasing=T)[1:(k+2)], 4)
round( sort(abs(x_mcp),decreasing=T)[1:(k+2)], 4)
round( sort(abs(x_scad),decreasing=T)[1:(k+2)], 4)

#support indices
sort(which(abs(x)>0))
sort(order(abs(sol_feas[,1]),decreasing=T)[1:k])
sort(order(abs(sol_bqn1[,1]),decreasing=T)[1:k])
sort(order(abs(x_lasso),decreasing=T)[1:k])
sort(order(abs(x_mcp),decreasing=T)[1:k])
sort(order(abs(x_scad),decreasing=T)[1:k])


#mean squared error
round( sqrt( mean( (sol[,1] - x)^2) ),5 )
round( sqrt( mean( (x_lasso - x)^2) ),5 )
round( sqrt( mean( (x_mcp - x)^2) ),5 )
round( sqrt( mean( (x_scad - x)^2) ),5 )
