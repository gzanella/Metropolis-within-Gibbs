#import libraries
#library(vioplot)
library(mcmcse)
library(matrixStats)
#library(R2OpenBUGS)
#library(rstan)
#library(rootSolve)
#simulate data
simulate_logit <- function(J,m, mu, tau, p){
  Y <- list()
  X <- list()
  for(j in 1:J){
    beta <- mu+1/sqrt(tau)*rnorm(p)
    x <- matrix(c(rep(1, m), runif((p-1)*m, -5, 5)), byrow = F, ncol = p, nrow = m)
    theta <- x%*%beta
    probs <- exp(theta)/(1+exp(theta))
    
    X[[j]] <- x
    Y[[j]] <- rbinom(m, size = 1, prob = probs)
  }
  return(list(X, Y))
}
source("Functions_HierCovariates.R")

#Simulation J = 30 and m = 30
m = 30
alpha <- beta <- 1
P <- c(1:5) #number of covariates
mu0 <- runif(P[length(P)], -1, 1)
tau0 <- rep(0.5, P[length(P)])

number_groups <- c(30)
replica <- 50
Iats_barker <- list()
Iats_barker5 <- list()
Iats_barker20 <- list()
Iats_independent <- list()
Iats_rwm <- list()

for(k in 1:length(P)){
  p <- P[k]
  mu <- mu0[1:p]
  tau <- tau0[1:p]
  Iats_barker[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_barker5[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_barker20[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_independent[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_rwm[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  
  burnin <- 1000
  iterations <- 1000
  #replicates
  for(i in 1:replica){
    print(paste("Replica", i, "with dimensionality", p))
    J_tot <- number_groups[length(number_groups)]
    data_set <- simulate_logit(J_tot,m, mu, tau, p)
    X_tot <- data_set[[1]]
    Y_tot <- data_set[[2]]
    
    #number of groups
    for(j in 1:length(number_groups)){
      J <- number_groups[j]
      Y <- Y_tot[1:J]
      X <- X_tot[1:J]
      
      
      #BARKER with 1 step
      sigma <- 2.4/(p^(1/6))
      reps = 1
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "barker")
      Iats_barker[[k]][i,j] <- out_MH[[5]]
      #print("Barker done")
      
      #BARKER with 5 steps
      sigma <- 2.4/(p^(1/6))
      reps = 5
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "barker")
      Iats_barker5[[k]][i,j] <- out_MH[[5]]
      
      #BARKER with 20 steps
      sigma <- 2.4/(p^(1/6))
      reps = 20
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "barker")
      Iats_barker20[[k]][i,j] <- out_MH[[5]]
      
      #INDEPENDENT
      sigma <- 2.4/(p^(1/2))
      reps = 1
      if(p <= 3){
        out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "independent")
        Iats_independent[[k]][i,j] <- out_MH[[5]]
      }
      
      
      #RWM
      sigma <- 2.4/(p^(1/2))
      reps <- 1
      if(p <= 3){
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "rwm")
      Iats_rwm[[k]][i,j] <- out_MH[[5]]
      }
      #print("Rwm done")
      
      #print
      to_print <- c(median(Iats_barker[[k]][1:i,j]), median(Iats_barker5[[k]][1:i,j]),median(Iats_barker20[[k]][1:i,j]),median(Iats_independent[[k]][1:i,j], na.rm = T), median(Iats_rwm[[k]][1:i,j]))
      names(to_print) <- c("Barker", "Barker5","Barker20", "Independent", "Rwm")
      print(round(to_print, 1))
      #cat("\n")
    }
  }
}


perf_barker1 <- sapply(Iats_barker, median)
perf_barker5 <- sapply(Iats_barker5, median)
perf_barker20 <- sapply(Iats_barker20, median)
perf_independent <- sapply(Iats_independent, median, na.rm = T)
perf_rwm <- sapply(Iats_rwm, median)



P <- 1:5
#logxy
Max <- max(c(perf_barker1[1:5],perf_barker5[1:5],perf_barker20[1:5], perf_rwm[1:3], perf_independent[1:3]))
Min <- min(c(perf_barker1[1:5],perf_barker5[1:5],perf_barker20[1:5], perf_rwm[1:3], perf_independent[1:3]))
plot(P[1:5], perf_barker1[1:5], type = "b", log = "",lwd = 2, pch = 16, col = gray(0.45), ylim = c(Min, Max), xlab = "N. covariates per group", ylab = "Int. aut. times ", cex = 1.5, cex.lab = 1.8, cex.axis = 2)
points(P[1:5], perf_barker5[1:5], type = "b", lwd = 2, pch = 17, col = gray(0.3), ylim = c(Min, Max),cex = 1.5)
points(P[1:5], perf_barker20[1:5], type = "b", lwd = 2, pch = 19, col = gray(0), ylim = c(Min, Max),cex = 1.5)
points(P[1:3], perf_independent[1:3], type = "b", lwd = 2, pch = 15, ylim = c(Min, Max), col = gray(0.6), cex = 1.5)
points(P[1:3], perf_rwm[1:3], type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.75), cex = 1.5)

plot(P, perf_barker1, type = "b", log = "",lwd = 2, pch = 16, col = gray(1), ylim = c(Min, Max), xlab = "Number of covariates", ylab = "Int. aut. times ", cex = 1.5, cex.lab = 1.8, cex.axis = 2)
legend(1,175, legend = c("MwG (Barker, 20 steps)", "MwG (Barker, 5 steps)", "MwG (Barker, 1 step)", "MwG (RWM, 1 step)", "MwG (IMH)"), lwd = 2, pch = c(19,17,16,18,15),
       col = c(gray(0), gray(0.3), gray(0.45), gray(0.75), gray(0.6)), bty = "n", y.intersp = 0.2)
