#computing logposterior
lpost_nonC <- function(y,x, beta, tau){
  p <- length(beta)
  theta <- 1/sqrt(tau)*beta
  prod_xtheta <- x%*%theta
  
  res <- sum(y*prod_xtheta-log(1+exp(prod_xtheta))) #likelihood
  res <- res-1/2*p*sum(beta^2) #prior of theta
  res <- as.vector(res)
  
  return(res)
}
grad_lpost_nonC <- function(y,x,beta, tau){
  p <- length(beta)
  theta <- 1/sqrt(tau)*beta
  prod_xtheta <- x%*%theta

  res <- 1/sqrt(tau)*matrix(y-exp(prod_xtheta)/(1+exp(prod_xtheta)), nrow = 1)%*%x #likelihood
  res <- res - p*beta
  res <- as.vector(res)

  return(res)
}

#GLOBAL PARAMETERS
sample_global_tau_nonC <- function(y, x, beta, tau, alpha1, alpha2, sigma_tau){
  log_tau<-log(tau)
  #proposal
  proposal <- log_tau+sigma_tau*rnorm(1)
  #acceptance probability
  log_acc <- lpost_nonC(y, x, beta, exp(proposal)) - lpost_nonC(y, x, beta, exp(log_tau))
  log_acc <- log_acc + dgamma(exp(proposal), alpha1, alpha2, log = T)-dgamma(exp(log_tau), alpha1, alpha2, log = T)
  log_acc <- log_acc +proposal-log_tau
  if(log(runif(1)) < log_acc){
    tau <- exp(proposal)
  }

  return(tau)
}

sample_global_tau_update_nonC <- function(y, x, beta, tau, alpha1, alpha2, sigma_tau, t){
  
  log_tau<-log(tau)
  #proposal
  proposal <- log_tau+sigma_tau*rnorm(1)
  #acceptance probability
  log_acc <- lpost_nonC(y, x, beta, exp(proposal)) - lpost_nonC(y, x, beta, exp(log_tau))
  log_acc <- log_acc + dgamma(exp(proposal), alpha1, alpha2, log = T)-dgamma(exp(log_tau), alpha1, alpha2, log = T)
  log_acc <- log_acc +proposal-log_tau
  if(log(runif(1)) < log_acc){
    tau <- exp(proposal)
  }

  #proposal
  proposal <- tau+sigma_tau*rnorm(1)
  if(proposal <0){
    #update sigma
    acc <- 0
    sigma2_log <- log(sigma_tau^2)+(1+t)^(-0.7)*(acc-0.23)
    sigma_tau <- sqrt(exp(sigma2_log))
    
    return(list(tau, sigma_tau))
  }
  #acceptance probability
  log_acc <- lpost_nonC(y, x, beta, proposal) - lpost_nonC(y, x, beta, tau)
  log_acc <- log_acc + dgamma(proposal, alpha1, alpha2, log = T)-dgamma(tau, alpha1, alpha2, log = T)
  if(log(runif(1)) < log_acc){
    tau <- proposal
  }
  
  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma2_log <- log(sigma_tau^2)+(1+t)^(-0.7)*(acc-0.23)
  sigma_tau <- sqrt(exp(sigma2_log))
  
  return(list(tau, sigma_tau))
}

#BARKER
rbarker_diag_nonC<-function(beta,c,sigma,diag_sd){
  # z<-sigma*rnorm(n=length(c),mean = diag_sd,sd = 0.1*diag_sd)
  z<-sigma*rnorm(n=length(c),mean = 0,sd = diag_sd)
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  
  return(beta+z*b)
}
log_q_ratio_diag_barker_nonC<-function(x,y,grad_x, grad_y){
  beta1<-  c(-grad_y*(x-y))
  beta2<-  c(-grad_x*(y-x))
  
  if(any(is.na(grad_y))){
    return(-Inf)
  }
  
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
sample_groups_vect_barker_nonC <- function(beta, y, x, tau, sigma,theta_mean, diag_sd){
  #proposal
  grad_beta <- grad_lpost_nonC(y,x,beta, tau)
  
  barker <- rbarker_diag_nonC(beta,grad_beta,sigma, diag_sd)
  grad_barker <- grad_lpost_nonC(y,x,barker, tau)
  
  #acceptance probability
  log_acc <- lpost_nonC(y,x,barker, tau)-lpost_nonC(y,x,beta, tau)+
    log_q_ratio_diag_barker_nonC(beta,barker,grad_beta,grad_barker)
  
  if(log(runif(1)) <= log_acc){
    beta <- barker
  }
  
  return(beta)
}
sample_groups_vect_update_sigma_barker_nonC <- function(beta, y, x, tau, sigma,theta_mean, diag_sd, t){
  
  #proposal
  grad_beta <- grad_lpost_nonC(y,x,beta, tau)
  
  barker <- rbarker_diag_nonC(beta,grad_beta,sigma, diag_sd)
  grad_barker <- grad_lpost_nonC(y,x,barker, tau)
  
  #acceptance probability
  log_acc <- lpost_nonC(y,x,barker, tau)-lpost_nonC(y,x,beta, tau)+
    log_q_ratio_diag_barker_nonC(beta,barker,grad_beta,grad_barker)

  if(log(runif(1)) <= log_acc){
    beta <- barker
  }

  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma2_log <- log(sigma^2)+(1+t)^(-0.7)*(acc-0.4)
  sigma <- sqrt(exp(sigma2_log))
  
  theta_mean<-theta_mean+(1+t)^(-0.7)*(beta-theta_mean)
  
  diag_sd2<-diag_sd^2+(1+t)^(-0.7)*(c((beta-theta_mean)^2)-diag_sd^2)
  diag_sd <- sqrt(diag_sd2)
  
  to_return <- list(beta, sigma, theta_mean, diag_sd)#matrix(c(theta, sigma), ncol = 2, byrow = F)
  return(to_return)
}

#RWM
sample_groups_vect_rwm_nonC <- function(beta, y, x,tau, sigma, theta_mean, diag_sd){
  #proposal
  p <- length(beta)
  proposal <- beta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost_nonC(y,x,proposal, tau)-lpost_nonC(y,x,beta, tau)
  
  if(log(runif(1)) <= log_acc){
    beta <- proposal
  }
  
  return(beta)
}
sample_groups_vect_update_sigma_rwm_nonC <- function(beta, y,x, tau, sigma, theta_mean, diag_sd, t){
  #proposal
  p <- length(beta)
  proposal <- beta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost_nonC(y,x,proposal, tau)-lpost_nonC(y,x,beta, tau)
  
  if(log(runif(1)) <= log_acc){
    beta <- proposal
  }
  
  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma2_log <- log(sigma^2)+(1+t)^(-0.7)*(acc-0.23)
  sigma <- sqrt(exp(sigma2_log))
  
  theta_mean<-theta_mean+(1+t)^(-0.7)*(beta-theta_mean)
  
  diag_sd2<-diag_sd^2+(1+t)^(-0.7)*(c((beta-theta_mean)^2)-diag_sd^2)
  diag_sd <- sqrt(diag_sd2)
  
  to_return <- list(beta, sigma, theta_mean, diag_sd)#matrix(c(theta, sigma), ncol = 2, byrow = F)
  return(to_return)
}


#METROPOLIS HASTINGS VECTORIZED
MH_sampler_vect_nonC <- function(Y,X, iterations, burnin, alpha1, alpha2, sigma = 2.4, tau0 = 1, adaptive = F, reps = 1, verbose = F, method = "barker"){
  #Y = vector of observatios
  #X = matrix of covariates
  #alpha1, alpha2 = prior parameters for tau
  #sigma = stepsize barker/rwm
  #tau0 = true parameter
  #reps = number of iterations MH per group
  #method = barker/rwm
  
  n <- length(Y)
  p <- ncol(X) #number of variables
  
  #starting point
  tau <- tau0 + runif(1, min = -1/sqrt(p), max = 1/sqrt(p))
  it_warm <- 10^4 #number iterations warm start
  
  #preconditioning
  sigma <- sigma
  sigma_tau <- 2.4
  theta_mean <- rep(0, p)
  diag_sd <- rep(1, p)
  
  beta <- 1/sqrt(p)*rnorm(p) #initialization
  sigma_store <- rep(0, it_warm)
  sigma_tau_store <- rep(0, it_warm)
  
  #BARKER
  if(method == "barker"){
    #warm start
    for(it in 1:(it_warm)){
      aux <- sample_groups_vect_update_sigma_barker_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd, it)
      beta <- aux[[1]]
      sigma <- aux[[2]]
      sigma_store[it] <- sigma
      theta_mean <- aux[[3]]
      diag_sd <- aux[[4]]
    }
    for(it in 1:(it_warm)){
      aux <- sample_global_tau_update_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau, it)
      tau <- aux[[1]]
      sigma_tau <- aux[[2]]
      sigma_tau_store[it] <- sigma_tau
    }
    if(verbose){
      print(paste("After warm:",tau))
      print(sigma)
      print(sigma_tau)
      plot(sigma_store, type = "l", main = "Sigma")
      plot(sigma_tau_store, type = "l", main = "Sigma Tau")
      cat("\n")
    }
    
    #burnin
    for(it in 1:burnin){
      for(k in 1:reps){
        #sample group parameters
        beta <- sample_groups_vect_barker_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau)
    }
    #iterations
    Tau <- rep(0, iterations)
    Beta <- matrix(0, nrow = iterations, ncol = p)
    for(it in 1:iterations){
      for(k in 1:reps){
        #sample group parameters
        beta <- sample_groups_vect_barker_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau)
      
      #save values
      Tau[it] <- tau
      Beta[it,] <- beta
    }
  }
  
  #RWM
  if(method == "rwm"){
    #warm start
    for(it in 1:(it_warm)){
      aux <- sample_groups_vect_update_sigma_rwm_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd, it)
      beta <- aux[[1]]
      sigma <- aux[[2]]
      theta_mean <- aux[[3]]
      diag_sd <- aux[[4]]
    }
    for(it in 1:(it_warm)){
      aux <- sample_global_tau_update_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau, it)
      tau <- aux[[1]]
      sigma_tau <- aux[[2]]
    }
    if(verbose){
      print(paste("After warm:",tau))
      print(sigma)
      print(sigma_tau)
      cat("\n")
    }
    
    #burnin
    for(it in 1:burnin){
      for(k in 1:reps){
        #sample group parameters
        beta <- sample_groups_vect_rwm_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau)
    }
    #iterations
    Tau <- rep(0, iterations)
    Beta <- matrix(0, nrow = iterations, ncol = p)
    for(it in 1:iterations){
      for(k in 1:reps){
        #sample group parameters
        beta <- sample_groups_vect_rwm_nonC(beta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau_nonC(Y, X, beta, tau, alpha1, alpha2, sigma_tau)
      
      #save values
      Tau[it] <- tau
      Beta[it,] <- beta
    }
  }
  
  ESS_min <- c(ess(Tau),ess(rowSums2(Beta^2)), min(ess(Beta)), min(ess(Beta^2)))
  IAT_max <- iterations/ESS_min
  if(verbose){
    print("Iats for different parameters:")
    print(c(IAT_max))
    cat("\n")
  }
  
  to_return <- list()
  to_return[[1]] <- Tau
  to_return[[2]] <- Beta
  to_return[[3]] <- sigma
  to_return[[4]] <- IAT_max
  
  return(to_return)
}