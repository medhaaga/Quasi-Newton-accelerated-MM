#################################################
######## Generalised Eigenvalues ################
#################################################

rm(list = ls())
library(pracma)
library(turboEM)
library(quasiNewtonMM)
library(daarem)
source("functions.R")

set.seed(1)
dim <- 500
C <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
D <- matrix(runif(dim^2, min = -5, max = 5), nrow = dim, ncol = dim)
A <- C + t(C)
B <- D %*% t(D)

N <- 10
start_rep <- matrix(rnorm(N*dim, mean = 0, sd = 10), nrow = N, ncol = dim)
tol <- 1e-3
D <- c("descent", "ascent")

for (d in 1:2)
{
  dir <- D[d]

  ###########################################
  ## Unaccelerated MM Algorithm
  ###########################################

  time_mm <- rep(NA, N)
  obj_mm <- rep(NA, N)
  eval_mm <- rep(NA, N)

   for (i in 1:N){
     print(i)
     now <- start_rep[i,]
     new <- start_rep[i,]
     diff <- 100

     iter <- 0
     start.time <- Sys.time()
     while((diff > tol))
     {
       iter <- iter + 1
       new <- update(now, A, B, dir = dir)
       diff <- sqrt(crossprod(new-now))
       now <- new
     }
     end.time <- Sys.time()

     time_mm[i] <- end.time - start.time
     obj_mm[i] <- rayleigh(new, A, B, dir)
     eval_mm[i] <- iter
   }
#
   print(quantile(time_mm, probs = c(.5, .25, .75)))
   print(quantile(eval_mm, probs = c(.5, .25, .75)))
   print(quantile((-1)^(d+1)*obj_mm, probs = c(.5, .25, .75)))

########################################
  ## BQN, q=1
  ########################################

  time_bqn1 <- rep(NA, N)
  obj_bqn1 <- rep(NA, N)
  eval_bqn1 <- rep(NA, N)

  for (i in 1:N){
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
              control = list(qn=1, tol = tol, maxiter = 5e4, objfn.inc=0.1))
    end.time <- Sys.time()

    if(fp$convergence){
      time_bqn1[i] <- end.time - start.time
      obj_bqn1[i] <- round(fp$value.objfn, 4)
      eval_bqn1[i] <- fp$fpevals
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_bqn1))))
  print(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)))
  print(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)))

#   ########################################
#   ## BQN, q=2
#   ########################################
# 
  time_bqn2 <- rep(NA, N)
  obj_bqn2 <- rep(NA, N)
  eval_bqn2 <- rep(NA, N)

  for (i in 1:N){
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
              control = list(qn=2, tol = tol, maxiter = 5e4, objfn.inc=0.001))
    end.time <- Sys.time()

    if(fp$convergence){
      time_bqn2[i] <- end.time - start.time
      obj_bqn2[i] <- fp$value.objfn
      eval_bqn2[i] <- fp$fpevals
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_bqn2))))
  print(quantile(time_bqn2[!is.na(time_bqn2)], probs = c(.5, .25, .75)))
  print(quantile(eval_bqn2[!is.na(eval_bqn2)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_bqn2[!is.na(obj_bqn2)], probs = c(.5, .25, .75)))

#   ########################################
#   ## BQN, q=5
#   ########################################

  time_bqn3 <- rep(NA, N)
  obj_bqn3 <- rep(NA, N)
  eval_bqn3 <- rep(NA, N)

  for (i in 1:N){
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- BQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
              control = list(qn=5, tol = tol, maxiter = 5e4, objfn.inc=1))
    end.time <- Sys.time()

    if(fp$convergence){
      time_bqn3[i] <- end.time - start.time
      obj_bqn3[i] <- fp$value.objfn
      eval_bqn3[i] <- fp$fpevals
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_bqn3))))
  print(quantile(time_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)))
  print(quantile(eval_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)))

  #######################################
  ## L-BQN
  ########################################

  print("Running LBQN")
  time_lbqn <- rep(NA, N)
  obj_lbqn <- rep(NA, N)
  eval_lbqn <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- LBQN(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
               control = list(m=10, tol = tol, maxiter = 5e4, objfn.inc=0.01))
    end.time <- Sys.time()

    if(fp$convergence){
      time_lbqn[i] <- end.time - start.time
      obj_lbqn[i] <- fp$value.objfn
      eval_lbqn[i] <- fp$fpevals
    }

  }
  print(paste("Number of failures: ", sum(is.na(time_lbqn))))
  print(quantile(time_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)))
  print(quantile(eval_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)))

  ##########################################
  ### SQUAREM v1
  #############################################

  print("Running SQUAREM-1")
  
  time_sq1 <- rep(NA, N)
  obj_sq1 <- rep(NA, N)
  eval_sq1 <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "squarem", control.method = list(K=1, version=1, objfn.inc=0.01),
                  control.run = list(tol = tol, maxiter = 3e4))
    end.time <- Sys.time()

    if(fp$convergence){
      time_sq1[i] <- end.time - start.time
      obj_sq1[i] <- fp$value.objfn
      eval_sq1[i] <- fp$fpeval
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_sq1))))
  print(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)))
  print(quantile(eval_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)))

  ##########################################
  ### SQUAREM v2
  #############################################

  print("Running SQUAREM-2")
  time_sq2 <- rep(NA, N)
  obj_sq2 <- rep(NA, N)
  eval_sq2 <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "squarem", control.method = list(K=1, version=2, objfn.inc=0.01),
                  control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()

    if(fp$convergence){
      time_sq2[i] <- end.time - start.time
      obj_sq2[i] <- fp$value.objfn
      eval_sq2[i] <- fp$fpeval
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_sq2))))
  print(quantile(time_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)))
  print(quantile(eval_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)))


  ##########################################
  ### SQUAREM v3
  #############################################

  print("Running SQUAREM-3")
  time_sq3 <- rep(NA, N)
  obj_sq3 <- rep(NA, N)
  eval_sq3 <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "squarem", control.method = list(K=1, version=3, objfn.inc=0.01),
                  control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()

    if(fp$convergence){
      time_sq3[i] <- end.time - start.time
      obj_sq3[i] <- fp$value.objfn
      eval_sq3[i] <- fp$fpeval
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_sq3))))
  print(quantile(time_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)))
  print(quantile(eval_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)))
  print(quantile((-1)^(d+1)*obj_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)))

  ##########################################
  ## ZAL for q=1
  ##########################################

  
  print("Running ZAL-1")
  time_zal <- rep(NA, N)
  obj_zal <- rep(NA, N)
  eval_zal <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "qn", control.method = list(qn=1),
                  control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    if(fp$convergence){
      time_zal[i] <- end.time - start.time
      obj_zal[i] <- fp$value.objfn
      eval_zal[i] <- fp$fpeval
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_zal))))
  print(quantile(time_zal[!is.na(time_zal)]))
  print(quantile(eval_zal[!is.na(time_zal)]))
  print(quantile((-1)^(d+1)*obj_zal[!is.na(time_zal)]))

  ##########################################
  ## ZAL for q=2
  ##########################################

  print("Running ZAL-2")
  time_zal2 <- rep(NA, N)
  obj_zal2 <- rep(NA, N)
  eval_zal2 <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "qn", control.method = list(qn=2),
                  control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()
    if(fp$convergence){
      time_zal2[i] <- end.time - start.time
      obj_zal2[i] <- fp$value.objfn
      eval_zal2[i] <- fp$fpeval
    }
  }
  print(sum(is.na(time_zal2)))
  print(paste("Number of failures: ", sum(is.na(time_zal2))))
  print(quantile(time_zal2[!is.na(time_zal2)]))
  print(quantile(eval_zal2[!is.na(time_zal2)]))
  print(quantile((-1)^(d+1)*obj_zal2[!is.na(time_zal2)]))

  ##########################################
  ## ZAL for q=min(p,10)
  ##########################################

  print("Running ZAL-10")
  time_zal3 <- rep(NA, N)
  obj_zal3 <- rep(NA, N)
  eval_zal3 <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- turboem(par = start, fixptfn = update, objfn = rayleigh, A=A, B=B, dir=dir,
                  method = "qn", control.method = list(qn=5),
                  control.run = list(tol = tol, maxiter = 5e4))
    end.time <- Sys.time()

    if(fp$convergence){
      time_zal3[i] <- end.time - start.time
      obj_zal3[i] <- fp$value.objfn
      eval_zal3[i] <- fp$fpeval
    }
  }

  print(paste("Number of failures: ", sum(is.na(time_zal3))))
  print(quantile(time_zal3))
  print(quantile(eval_zal3))
  print(quantile((-1)^(d+1)*obj_zal3))

  ##########################################
  ## DAAREM
  ##########################################

  print("Running DAAREM")
  time_dar <- rep(NA, N)
  obj_dar <- rep(NA, N)
  eval_dar <- rep(NA, N)

  for (i in 1:N){
    print(i)
    start <- start_rep[i,]
    start.time <- Sys.time()
    fp <- daarem(par = start, fixptfn = update, objfn = daarem.objective, A=A, B=B,
                 dir=dir, control = list(tol = tol, maxiter = 5e4, mon.tol=0.01))
    end.time <- Sys.time()

    time_dar[i] <- end.time - start.time
    obj_dar[i] <- -fp$value.objfn
    eval_dar[i] <- fp$fpeval

  }

  print(quantile(round(time_dar[1], 3)))
  print(quantile(round(eval_dar[1], 3)))
  print(quantile(round((-1)^(d+1)*obj_dar[1], 3)))

  save(time_qn1, time_bqn1, time_bqn2, time_bqn3, time_lbqn, time_sq1, time_sq2, time_sq3, time_zal, time_zal2, time_zal3, time_dar,
       eval_qn1, eval_bqn1, eval_bqn2, eval_bqn3, eval_lbqn, eval_sq1, eval_sq2, eval_sq3, eval_zal, eval_zal2, eval_zal3, eval_dar,
       obj_qn1, obj_bqn1, obj_bqn2, obj_bqn3, obj_lbqn, obj_sq1, obj_sq2, obj_sq3, obj_zal, obj_zal2, obj_zal3, obj_dar, file = paste("Out/eigen-objects_", dir, "_dim", dim, "_sd1e2.Rdata", sep = ""))

}


for (d in 1:2){
  dir = D[d]
  load(file = paste("Out/eigen-objects_", dir, "_sd1e2.Rdata", sep = ""))

  print(paste("Case = ", dir))
  print("MM")
  print(paste("Number of failures: ", sum(is.na(time_mm))))
  print(round(quantile(time_mm, probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_mm, probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_mm, probs = c(.5, .25, .75)),3))
  
  print("BQN-1")
  print(paste("Number of failures: ", sum(is.na(time_bqn1))))
  print(round(quantile(time_bqn1[!is.na(time_bqn1)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_bqn1[!is.na(eval_bqn1)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_bqn1[!is.na(obj_bqn1)], probs = c(.5, .25, .75)),3))
  
  print("BQN-2")
  print(paste("Number of failures: ", sum(is.na(time_bqn2))))
  print(round(quantile(time_bqn2[!is.na(time_bqn2)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_bqn2[!is.na(eval_bqn2)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_bqn2[!is.na(obj_bqn2)], probs = c(.5, .25, .75)),3))
  
  print("BQN-3")
  print(paste("Number of failures: ", sum(is.na(time_bqn3))))
  print(round(quantile(time_bqn3[!is.na(time_bqn3)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_bqn3[!is.na(eval_bqn3)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_bqn3[!is.na(obj_bqn3)], probs = c(.5, .25, .75)),3))
  
  print("LBQN")
  print(paste("Number of failures: ", sum(is.na(time_lbqn))))
  print(round(quantile(time_lbqn[!is.na(time_lbqn)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_lbqn[!is.na(eval_lbqn)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_lbqn[!is.na(obj_lbqn)], probs = c(.5, .25, .75)),3))
  
  print("SQUAREM-1")
  print(paste("Number of failures: ", sum(is.na(time_sq1))))
  print(round(quantile(time_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_sq1[!is.na(time_sq1)], probs = c(.5, .25, .75)),3))
  
  print("SQUAREM-2")
  print(paste("Number of failures: ", sum(is.na(time_sq2))))
  print(round(quantile(time_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_sq2[!is.na(time_sq2)], probs = c(.5, .25, .75)),3))
  
  print("SQUAREM-3")
  print(paste("Number of failures: ", sum(is.na(time_sq3))))
  print(round(quantile(time_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)),3))
  print(round(quantile(eval_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)),3))
  print(round(quantile((-1)^(d+1)*obj_sq3[!is.na(time_sq3)], probs = c(.5, .25, .75)),3))
  
  print("ZAL")
  print(paste("Number of failures: ", sum(is.na(time_zal))))
  print(round(quantile(time_zal[!is.na(time_zal)], probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile(eval_zal[!is.na(time_zal)], probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile((-1)^(d+1)*obj_zal[!is.na(time_zal)], probs = c(0.5, 0.25, 0.75)),3))
  
  print("ZAL-2")
  print(sum(is.na(time_zal2)))
  print(paste("Number of failures: ", sum(is.na(time_zal2))))
  print(round(quantile(time_zal2[!is.na(time_zal2)], probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile(eval_zal2[!is.na(time_zal2)], probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile((-1)^(d+1)*obj_zal2[!is.na(time_zal2)], probs = c(0.5, 0.25, 0.75)),3))
  
  print("ZAL-3")
  print(paste("Number of failures: ", sum(is.na(time_zal3))))
  print(round(quantile(time_zal3, probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile(eval_zal3, probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile((-1)^(d+1)*obj_zal3, probs = c(0.5, 0.25, 0.75)),3))
  
  print("DAAREM")
  print(paste("Number of failures: ", sum(is.na(time_dar))))
  print(round(quantile(time_dar, probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile(eval_dar, probs = c(0.5, 0.25, 0.75)),3))
  print(round(quantile(obj_dar, probs = c(0.5, 0.25, 0.75)),3))
}


