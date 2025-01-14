#import libraries

library(mcmcse)
library(matrixStats)

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


#Simulation p = 5 and m = 30
m = 30
alpha <- beta <- 1
P <- c(5) #number of covariates
mu0 <- runif(10, -1, 1) #true mean hyperparameter
tau0 <- c(2,1,1,3,2) #true precision hyperparameters

number_groups <- c(120, 160, 200, 240, 280)
replica <- 50
Iats_barker <- list()
Iats_barker5 <- list()
Iats_barker10 <- list()
Iats_rwm <- list()
#reps <- 1
for(k in 1:length(P)){
  p <- P[k]
  mu <- mu0[1:p]
  tau <- tau0[1:p]
  Iats_barker[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_barker5[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_barker10[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_rwm[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  
  burnin <- 1000
  iterations <- 1500
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
      
      
      #BARKER
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
      
      #BARKER with 10 steps
      sigma <- 2.4/(p^(1/6))
      reps = 10
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "barker")
      Iats_barker10[[k]][i,j] <- out_MH[[5]]
      
      #RWM
      sigma <- 2.4/(p^(1/2))
      reps <- 1
      out_MH <- MH_sampler_vect(Y, X, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = T, verbose = F, reps = reps, method = "rwm")
      Iats_rwm[[k]][i,j] <- out_MH[[5]]
      #print("Rwm done")
      
      #print
      to_print <- c(median(Iats_barker[[k]][1:i,j]), median(Iats_barker5[[k]][1:i,j]),median(Iats_barker10[[k]][1:i,j]), median(Iats_rwm[[k]][1:i,j]))
      names(to_print) <- c("Barker", "Barker5","Barker10", "Rwm")
      print(round(to_print, 1))
      #cat("\n")
    }
  }
}
k <- 1

perf_barker <- Iats_barker[[k]]
perf_barker5 <- Iats_barker5[[k]]
perf_barker10 <- Iats_barker10[[k]]
perf_rwm <- Iats_rwm[[k]]


#saveRDS(perf_barker, file = "barker_Hier5Covariates.rds")
#saveRDS(perf_barker5, file = "barker5_Hier5Covariates.rds")
#saveRDS(perf_barker10, file = "barker5_Hier5Covariates.rds")
#saveRDS(perf_rwm, file = "rwm_Hier5Covariates.rds")
#saveRDS(number_groups, file = "number_groups_Hier5Covariates.rds")
Med_barker <- apply(perf_barker, 2, median)
Med_barker5 <- apply(perf_barker5, 2, median)
Med_barker10 <- apply(perf_barker10, 2, median)
Med_rwm <- apply(perf_rwm, 2, median)


#logxy
Max <- max(c(Med_barker, Med_barker5,Med_barker10, Med_rwm))
Min <- min(c(Med_barker, Med_barker5,Med_barker10, Med_rwm))
plot(number_groups, Med_barker, type = "b", lwd = 2, pch = 19, col = gray(0.6), ylim = c(Min, Max),log = "xy", xlab = "Groups", ylab = "Int. aut. times ", cex = 1.5, cex.lab = 1.8, cex.axis = 2)
points(number_groups, Med_barker5, type = "b", lwd = 2, pch = 17, ylim = c(Min, Max), col = gray(0.3), cex = 1.5)
points(number_groups, Med_barker10, type = "b", lwd = 2, pch = 15, ylim = c(Min, Max), col = gray(0), cex = 1.5)
points(number_groups, Med_rwm, type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.8), cex = 1.5)

plot(number_groups, Med_barker, type = "b", lwd = 2, pch = 19, col = gray(1), ylim = c(Min, Max), xlab = "Groups", ylab = "Int. aut. times ", cex = 1.5, cex.lab = 1.8, cex.axis = 2)
legend(150,250, legend = c("MwG (Barker, 10 steps)", "MwG (Barker, 5 steps)", "MwG (Barker, 1 step)", "MwG (RWM, 1 step)"), lwd = 2, pch = c(19,17,15,18),
       col = c(gray(0), gray(0.3), gray(0.6), gray(0.8)), bty = "n", y.intersp = 0.2)
