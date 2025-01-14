library(mcmcse)
library(matrixStats)

#### SET TO TRUE TO DO A SHORT RUN WITH MORE MONTE CARLO ERROR; FALSE TO DO A LONG RUN WITH LESS MONTE CARLO ERROR
SHORT_RUN<-TRUE

#### SET WORKING DIRECTORY AND LOAD AUXILIARY FUNCTIONS
setwd("C:/Users/ZanellaG/Dropbox/ResearchProjects/Fil/Filippo/MwG/Code/Code_Figures_4_5_selected/")
source("Functions_LogisticRegression.R")
source("Functions_LogisticRegression_nonC.R")

### SET THE NUMBER OF DATAPOINTS AND PARAMETERS FOR THE SIMULATION
n_vec <- c(10,30,50)
p_vec <- 2*n_vec

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
    
    #Gibbs
    reps <- 100
    sigma <- 2.4/p^(1/6)
    out_MH <- MH_sampler_vect(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "barker")
    Iats_gibbs[i, k] <- max(out_MH[[4]])
    
    # non-centered parametrization
    
    #rwm
    sigma <- 2.4/p^(1/2)
    reps <- 1
    out_MH <- MH_sampler_vect_nonC(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "rwm")
    Iats_rwm_nonC[i, k] <- max(out_MH[[4]])
    
    #barker
    reps <- 1
    sigma <- 2.4/p^(1/6)
    out_MH <- MH_sampler_vect_nonC(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "barker")
    Iats_barker_nonC[i, k] <- max(out_MH[[4]])
    
    #Gibbs
    reps <- 100
    sigma <- 2.4/p^(1/6)
    out_MH <- MH_sampler_vect_nonC(Y,X, iterations, burnin, alpha1, alpha2, sigma = sigma, tau0 = tau, adaptive = F, reps = reps, verbose = F, method = "barker")
    Iats_gibbs_nonC[i, k] <- max(out_MH[[4]])
    
    print("Centered")
    to_print <- c(median(Iats_gibbs[1:i, k]))
    names(to_print) <- c("Gibbs")
    print(round(to_print, 1))
    print("Non-centered")
    to_print <- c(median(Iats_rwm_nonC[1:i, k]), median(Iats_barker_nonC[1:i, k]), median(Iats_gibbs_nonC[1:i, k]))
    names(to_print) <- c("Rwm", "Barker", "Gibbs")
    print(round(to_print, 1))
    cat("\n")
  }
}

### COMPUTE MEDIAN IATs
Iats_list <- vector("list", 6)
Iats_list[[3]] <- Iats_gibbs
Iats_list[[4]] <- Iats_rwm_nonC
Iats_list[[5]] <- Iats_barker_nonC
Iats_list[[6]] <- Iats_gibbs_nonC

Med_gibbs <- apply(Iats_list[[3]],2,median)
Med_rwm_nonC <- apply(Iats_list[[4]],2,median)
Med_barker_nonC <- apply(Iats_list[[5]],2,median)
Med_gibbs_nonC <- apply(Iats_list[[6]],2,median)

### PLOT FIGURE 5
Max <- max(Med_gibbs, Med_rwm_nonC, Med_barker_nonC, Med_gibbs_nonC)
Min <- min(Med_gibbs, Med_rwm_nonC, Med_barker_nonC, Med_gibbs_nonC)
plot(p_vec, Med_gibbs, type = "b",lty = 1, lwd = 2, pch = 19, 
     col = gray(0), ylim = c(Min, Max),
     log = "", xlab = "Number of covariates", ylab = "Int. aut. times ", 
     cex = 1.5, cex.lab = 1.8, cex.axis = 2)
points(p_vec, Med_rwm_nonC, type = "b",lty = 2, lwd = 2, pch = 17, col = gray(0.8), ylim = c(Min, Max),cex = 1.5)
points(p_vec, Med_barker_nonC, type = "b",lty = 2, lwd = 2, pch = 18, col = gray(0.4), ylim = c(Min, Max),cex = 1.5)
points(p_vec, Med_gibbs_nonC, type = "b",lty = 2, lwd = 2, pch = 19, col = gray(0), ylim = c(Min, Max),cex = 1.5)
legend(10,450, 
       legend = c("MwG (Barker, 100 steps) cent.", "MwG (Barker, 100 steps) non-cent.", 
                          "MwG (Barker, 1 step) non-cent.", "MwG (RWM, 1 step) non-cent."), 
       lwd = 2, pch = c(19,19,18, 17),
       col = c(gray(0), gray(0), gray(0.4), gray(0.8)), bty = "n", 
       y.intersp = 0.25, cex = 1.2, lty = c(1,2,2,2))