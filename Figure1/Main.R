#import libraries
library(mcmcse)
library(matrixStats)
library(rstan)



source("https://raw.githubusercontent.com/gzanella/Metropolis-within-Gibbs/refs/heads/main/Figure1/Functions_MH.R")
source("https://raw.githubusercontent.com/gzanella/Metropolis-within-Gibbs/refs/heads/main/Figure1/Functions_MALA.R")


##if you also want exact Gibbs implementation:
# library(ars)
# source("https://raw.githubusercontent.com/gzanella/Metropolis-within-Gibbs/refs/heads/main/Figure1/Functions_Gibbs.R")
## in the simulation Gibbs is replaced by MwG with large number of steps per iteration

#simulate data
simulate_logit <- function(J,m, mu, tau){
  #simulate data according to hierarchical logit model
  theta <- mu+1/sqrt(tau)*rnorm(J)
  probs <- exp(theta)/(1+exp(theta))
  Y <- rbinom(J, size = m, prob = probs)
  return(Y)
}



#FINAL SIMULATION
M <- c(10) #number of data per group
mu <- 1 #true mean
tau <- 1 #true precision
alpha <- beta <- 1 #hyperparameter Gamma prior on tau
burnin <- 2000
iterations <- 4000

number_groups <- 1*2^c(7:11)
replica <- 50 #replications of the experiment
#to store the IATs
Iats_Gibbs <- list()
Iats_MH <- list()
Iats_MALA <- list()
Iats_STAN <- list()
leaps_STAN_store <- list()
#parameters MALA
gamma0 <- 2.4
reps_MALA <- 1
sigma <- 1
#parameter Gibbs and MH
reps_gibbs <- 5
reps_MH <- 1
for(k in 1:length(M)){
  m <- M[k] #number of observations per group
  Iats_Gibbs[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_MH[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_MALA[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  Iats_STAN[[k]] <- matrix(0,nrow = replica, ncol = length(number_groups))
  leaps_STAN_store[[k]] <- matrix(0, ncol = length(number_groups), nrow = replica)
  for(i in 1:replica){
    #replicas
    print(paste("Replica", i, "with group size", m))
    J_tot <- number_groups[length(number_groups)]
    Y_tot <- simulate_logit(J_tot,m, mu, tau)
    for(j in 1:length(number_groups)){
      J <- number_groups[j]
      y <- Y_tot[1:J]
      
      
      #Gibbs: Mwg with mutiple steps per iteration
      out_gibbs <- MH_sampler_vect(y, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = F, reps = reps_gibbs)
      sim <- matrix(c(out_gibbs[[1]], log(out_gibbs[[2]]), out_gibbs[[3]]), nrow = iterations, byrow = F)
      sd_gibbs <- colSds(sim) #standard deviations, used to tune MALA
      Iats_Gibbs[[k]][i,j] <- out_gibbs[[5]]
      
        # Exact Gibbs
      #out_gibbs <- Gibbs_sampler(y, m, iterations, burnin, alpha, beta, mu0 = mu, tau0 = tau)
      #sim <- matrix(c(out_gibbs[[1]], log(out_gibbs[[2]]), out_gibbs[[3]]), nrow = iterations, byrow = F)
      #sd_gibbs <- colSds(sim) #standard deviations, used to tune MALA
      #Iats_Gibbs[[k]][i,j] <- out_gibbs[[4]]
      
      
      
      #MH
      out_MH <- MH_sampler_vect(y, m, iterations, burnin, alpha, beta, sigma, mu0 = mu, tau0 = tau, adaptive = F, reps = reps_MH)
      Iats_MH[[k]][i,j] <- out_MH[[5]]
      #print("MwG done")
      
      #MALA
      out_MALA <- MALA_sampler(y, m, iterations, burnin, alpha, beta, gamma0 = gamma0,sd = sd_gibbs, mu0 = mu, tau0 = tau, reps = reps_MALA)
      Iats_MALA[[k]][i,j] <- out_MALA[[3]]
      #print("MALA done")
      
      
      
      #STAN
      sink(file = "a", type = c("output", "message"))
      out_STAN <- suppressMessages(suppressWarnings(stan(file="hierarchical_logit_stan.stan",algorithm = "NUTS", refresh = -1,data=list(m = m,J = J,y = y,alpha = alpha, beta = beta), chains = 1, iter = iterations+burnin, warmup = burnin))) 
      sink(type="output")
      #print("Stan done")
      
      #compute number of leaps per iteration
      draws_stan <- extract(out_STAN, permuted = F)[,1,]
      draws_stan <- draws_stan[,-ncol(draws_stan)]
      perf_STAN <- leap_ESS(out_STAN,"NUTS",t_burn = burnin,print=F)
      ess_STAN <- perf_STAN[[2]][-length(perf_STAN[[2]])]
      leaps_STAN <- perf_STAN[[1]]
      Iats_STAN[[k]][i,j] <- iterations/min(ess_STAN)
      leaps_STAN_store[[k]][i,j] <- leaps_STAN
      
      
      #print
      to_print <- c(Iats_Gibbs[[k]][i,j], Iats_MH[[k]][i,j], Iats_MALA[[k]][i,j], Iats_STAN[[k]][i,j]*leaps_STAN)
      names(to_print) <- c("Gibbs", "MwG", "MALA", "Stan")
      print(round(to_print, 2))
      #cat("\n")
    }
  }
}
#summaries
perf_gibbs <- Iats_Gibbs[[1]]
perf_MH <- Iats_MH[[1]]
perf_MALA <- Iats_MALA[[1]]
perf_STAN <- Iats_STAN[[1]]
leaps_STAN <- leaps_STAN_store[[1]]

#saveRDS(perf_gibbs, file = "gibbs10.rds")
#saveRDS(perf_MH, file = "MH10.rds")
#saveRDS(perf_MALA, file = "MALA10.rds")
#saveRDS(perf_STAN, file = "STAN10.rds")
#saveRDS(leaps_STAN, file = "leaps10.rds")
#saveRDS(number_groups, file = "number_groups10.rds")


Med_gibbs <- apply(perf_gibbs, 2, median)
Med_MH <- apply(perf_MH, 2, median)
Med_MALA <- apply(perf_MALA, 2, median)
Med_STAN <- apply(perf_STAN, 2, median)
Med_leaps <- apply(leaps_STAN, 2, median)
Med_STAN_leaps <- apply(perf_STAN*leaps_STAN, 2, median)




#plots
#logxy
Max <- max(c(Med_gibbs, Med_MH, Med_MALA, Med_STAN, Med_STAN_leaps))
Min <- min(c(Med_gibbs, Med_MH, Med_MALA, Med_STAN, Med_STAN_leaps))
plot(number_groups, Med_gibbs, type = "b", lwd = 2, pch = 19, col = gray(0), xlim = c(number_groups[1], number_groups[length(number_groups)]), ylim = c(Min, Max),log = "xy", xlab = "Groups", ylab = "IAT times lik. evaluations", cex = 1.5, cex.lab = 1.4, cex.axis = 2)
points(number_groups, Med_MH, type = "b", lwd = 2, pch = 17, ylim = c(Min, Max), col = gray(0.3), cex = 1.5)
points(number_groups, Med_MALA, type = "b", lwd = 2, pch = 15, ylim = c(Min, Max), col = gray(0.5), cex = 1.5)
points(number_groups, Med_STAN_leaps, type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.7), cex = 1.5)

#legend
plot(1:10, 1:10, type = "b", lwd = 2, pch = 15, col = "white", cex = 1.5, xlab = "", ylab = "")
legend(2,8, legend = c("Gibbs", "Metropolis-within-Gibbs", "Optimally tuned MALA", "Stan HMC"), lwd = 2, pch = c(19,17,15,18),
       col = c(gray(0), gray(0.3), gray(0.5), gray(0.7)), bty = "n", y.intersp = 2)
#y.intersp = 0.0001
