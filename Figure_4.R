library(mcmcse)
library(matrixStats)

#### SET TO TRUE TO DO A SHORT RUN WITH MORE MONTE CARLO ERROR; FALSE TO DO A LONG RUN WITH LESS MONTE CARLO ERROR
SHORT_RUN<-TRUE

#### SET WORKING DIRECTORY AND LOAD AUXILIARY FUNCTIONS
setwd("C:/Users/ZanellaG/Dropbox/ResearchProjects/Fil/Filippo/MwG/Code/Code_Figures_4_5_selected/")
source("Functions_LogisticRegression.R")
source("Functions_LogisticRegression_nonC.R")

### SET THE NUMBER OF DATAPOINTS AND PARAMETERS FOR THE SIMULATION
n_vec <- c(20,40,60,80,100)
p_vec <- rep(5, length(n_vec))

#### SET LENGTH OF CHAINS AND NUMBER OF MONTE CARLO REPLICAS
if(SHORT_RUN){
  burnin <- 10^2
  iterations <- 10^3
  replica <- 3
}else{
  burnin <- 500
  iterations <- 10^5
  replica <- 50
}

### SET THE PARAMETERS OF THE GAMMA(alpha1,alpha2) PRIOR FOR TAU
alpha1 <- alpha2 <- 1

### SET THE DATA-GENERATING TAU
tau <- 1

### INITIALIZE MATRICES
Iats_rwm <- matrix(0, nrow = replica, ncol = length(n_vec))
Iats_barker <- matrix(0, nrow = replica, ncol = length(n_vec))
Iats_gibbs <- matrix(0, nrow = replica, ncol = length(n_vec))

Iats_rwm_nonC <- matrix(0, nrow = replica, ncol = length(n_vec))
Iats_barker_nonC <- matrix(0, nrow = replica, ncol = length(n_vec))
Iats_gibbs_nonC <- matrix(0, nrow = replica, ncol = length(n_vec))

### RUN SIMULATIONS
for(k in 1:length(n_vec)){
  n <- n_vec[k]
  p <- p_vec[k]
  for(i in 1:replica){
    print(paste("Replica", i, "with size", n, "and dimensionality", p))
    data_set <- simulate_logit(p, n, tau)
    Y <- data_set[[2]]
    X <- data_set[[1]]
    
    # centered parametrization
    
    #rwm
    sigma <- 2.4/p^(1/2)
    reps <- 1
    out_MH <- MH_sampler_vect(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "rwm")
    Iats_rwm[i, k] <- max(out_MH[[4]])
    
    #barker
    reps <- 1
    sigma <- 2.4/p^(1/6)
    out_MH <- MH_sampler_vect(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "barker")
    Iats_barker[i, k] <- max(out_MH[[4]])
    
    #Gibbs
    reps <- 100
    sigma <- 2.4/p^(1/6)
    out_MH <- MH_sampler_vect(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "barker")
    Iats_gibbs[i, k] <- max(out_MH[[4]])
    
    print("Centered")
    to_print <- c(median(Iats_rwm[1:i, k]), median(Iats_barker[1:i, k]), median(Iats_gibbs[1:i, k]))
    names(to_print) <- c("Rwm", "Barker", "Gibbs")
    print(round(to_print, 1))
    cat("\n")
  }
}

### COMPUTE MEDIAN IATs
Iats_list <- vector("list", 6)
Iats_list[[1]] <- Iats_rwm
Iats_list[[2]] <- Iats_barker
Iats_list[[3]] <- Iats_gibbs

Med_rwm <- apply(Iats_list[[1]],2,median)
Med_barker <- apply(Iats_list[[2]],2,median)
Med_gibbs <- apply(Iats_list[[3]],2,median)

### PLOT FIGURE 4
Max <- max(Med_gibbs, Med_rwm, Med_barker)
Min <- min(Med_gibbs, Med_rwm, Med_barker)
plot(n_vec, Med_rwm, type = "b", lwd = 2, pch = 17, col = gray(0.8),
     ylim = c(Min, Max), xlab = "Number of observations", 
     ylab = "Int. aut. times ", cex = 1.5, cex.lab = 1.8, cex.axis = 2)
points(n_vec, Med_gibbs, type = "b",lty = 1, lwd = 2, pch = 19,col = gray(0), cex = 1.5)
points(n_vec, Med_barker, type = "b",lty = 1, lwd = 2, pch = 18, col = gray(0.4), cex = 1.5)
legend("right",#10,450,
       legend = c("MwG (Barker, 100 steps)",
                  "MwG (Barker, 1 steps)",
                  "MwG (RWM, 1 step)"), 
       lwd = 2, pch = c(19, 18, 17),
       col = c(gray(0), gray(0.4), gray(0.8)), bty = "n", 
       y.intersp = 1,#0.25, 
       cex = 1.2, lty = c(1,1,1))



