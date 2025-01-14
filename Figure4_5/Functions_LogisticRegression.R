simulate_logit <- function(p,n,tau){
  #covariates
  X <- matrix(c(rep(1, n), runif((p-1)*n, -5, 5)), byrow = F, ncol = p, nrow = n)
  
  #simulate alpha2
  alpha2 <- 1/sqrt(p*tau)*rnorm(p)
  
  #probabilities
  theta <- X%*%alpha2
  probs <- exp(theta)/(1+exp(theta))
  
  
  Y <- rbinom(n, size = 1, prob = probs)
  
  return(list(X, Y))
}

#computing logposterior
lpost <- function(y,x, theta, tau){
  p <- length(theta)
  prod_xtheta <- x%*%theta
  
  res <- sum(y*prod_xtheta-log(1+exp(prod_xtheta))) #likelihood
  res <- res-tau/2*p*sum(theta^2) #prior of theta
  res <- as.vector(res)
  
  return(res)
}
grad_lpost <- function(y,x,theta, tau){
  p <- length(theta)
  prod_xtheta <- x%*%theta

  res <- matrix(y-exp(prod_xtheta)/(1+exp(prod_xtheta)), nrow = 1)%*%x #likelihood
  #res <- t(res)-tau*Prec%*%theta #prior on theta
  res <- res - tau*p*theta
  res <- as.vector(res)

  return(res)
}

#GLOBAL PARAMETERS
sample_global_tau <- function(theta, alpha1, alpha2){
  p <- length(theta)
  res <- rgamma(1, alpha1+p/2, alpha2+0.5*p*sum(theta^2))
  
  return(res)
}

#BARKER
rbarker_diag<-function(theta,c,sigma,diag_sd){
  # z<-sigma*rnorm(n=length(c),mean = diag_sd,sd = 0.1*diag_sd)
  z<-sigma*rnorm(n=length(c),mean = 0,sd = diag_sd)
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  return(theta+z*b)
}
log_q_ratio_diag_barker<-function(x,y,grad_x, grad_y){
  if(any(is.na(grad_y))){
    return(-Inf)
  }
  alpha21<-  c(-grad_y*(x-y))
  alpha22<-  c(-grad_x*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(alpha21,0)+log1p(exp(-abs(alpha21))))+
      (pmax(alpha22,0)+log1p(exp(-abs(alpha22))))
  ))
}
sample_groups_vect_barker <- function(theta, y, x, tau, sigma,theta_mean, diag_sd){
  #proposal
  grad_theta <- grad_lpost(y,x,theta, tau)
  
  barker <- rbarker_diag(theta,grad_theta,sigma, diag_sd)
  grad_barker <- grad_lpost(y,x,barker, tau)
  
  #acceptance probability
  log_acc <- lpost(y,x,barker, tau)-lpost(y,x,theta, tau)+
    log_q_ratio_diag_barker(theta,barker,grad_theta,grad_barker)
  # print(log_q_ratio_diag_barker(theta,barker,grad_theta,grad_barker))
  
  if(log(runif(1)) <= log_acc){
    theta <- barker
    #print("Accepted")
  }
  
  return(theta)
}
sample_groups_vect_update_sigma_barker <- function(theta, y, x, tau, sigma,theta_mean, diag_sd, t){
  
  #proposal
  grad_theta <- grad_lpost(y,x,theta, tau)
  
  barker <- rbarker_diag(theta,grad_theta,sigma, diag_sd)
  grad_barker <- grad_lpost(y,x,barker, tau)
  
  #acceptance probability
  log_acc <- lpost(y,x,barker, tau)-lpost(y,x,theta, tau)+
    log_q_ratio_diag_barker(theta,barker,grad_theta,grad_barker)
  
  if(log(runif(1)) <= log_acc){
    theta <- barker
    #print("Accepted")
  }
  
  
  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma2_log <- log(sigma^2)+(1+t)^(-0.7)*(acc-0.4)
  sigma <- sqrt(exp(sigma2_log))
  
  theta_mean<-theta_mean+(1+t)^(-0.7)*(theta-theta_mean)
  
  diag_sd2<-diag_sd^2+(1+t)^(-0.7)*(c((theta-theta_mean)^2)-diag_sd^2)
  diag_sd <- sqrt(diag_sd2)
  
  #print(summary(sigma))
  to_return <- list(theta, sigma, theta_mean, diag_sd)#matrix(c(theta, sigma), ncol = 2, byrow = F)
  return(to_return)
}

#RWM
sample_groups_vect_rwm <- function(theta, y, x,tau, sigma, theta_mean, diag_sd){
  #proposal
  p <- length(theta)
  proposal <- theta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost(y,x,proposal, tau)-lpost(y,x,theta, tau)
  
  if(log(runif(1)) <= log_acc){
    theta <- proposal
    #print("Accepted")
  }
  
  return(theta)
}
sample_groups_vect_update_sigma_rwm <- function(theta, y,x, tau, sigma, theta_mean, diag_sd, t){
  #proposal
  p <- length(theta)
  proposal <- theta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost(y,x,proposal, tau)-lpost(y,x,theta, tau)
  
  if(log(runif(1)) <= log_acc){
    theta <- proposal
    #print("Accepted")
  }
  
  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma2_log <- log(sigma^2)+(1+t)^(-0.7)*(acc-0.23)
  sigma <- sqrt(exp(sigma2_log))
  
  theta_mean<-theta_mean+(1+t)^(-0.7)*(theta-theta_mean)
  
  diag_sd2<-diag_sd^2+(1+t)^(-0.7)*(c((theta-theta_mean)^2)-diag_sd^2)
  diag_sd <- sqrt(diag_sd2)
  
  #print(summary(sigma))
  to_return <- list(theta, sigma, theta_mean, diag_sd)#matrix(c(theta, sigma), ncol = 2, byrow = F)
  return(to_return)
}


#METROPOLIS HASTINGS VECTORIZED
MH_sampler_vect <- function(Y,X, iterations, burnin, alpha1, alpha2, sigma = 2.4, tau0 = 1, adaptive = F, reps = 1, verbose = F, method = "barker"){
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
  it_warm <- 10^4#as.integer(100*log(p))+1 #number iterations warm start
  
  #preconditioning
  sigma <- sigma
  theta_mean <- rep(0, p)
  diag_sd <- rep(1, p)
  
  theta <- 1/sqrt(p)*rnorm(p) #initialization
  sigma_store <- rep(0, it_warm)
  
  #BARKER
  if(method == "barker"){
    #warm start
    for(it in 1:(it_warm)){
      aux <- sample_groups_vect_update_sigma_barker(theta, Y, X, tau, sigma,theta_mean, diag_sd, it)
      theta <- aux[[1]]
      sigma <- aux[[2]]
      sigma_store[it] <- sigma
      theta_mean <- aux[[3]]
      diag_sd <- aux[[4]]
    }
    if(verbose){
      print("Adaptively tuned parameters:")
      print(sigma)
      print(c(min(diag_sd), max(diag_sd)))
      plot(sigma_store, type = "l", main = "Adaptive sigma")
      cat("\n")
    }
    
    #burnin
    
    #IF ADAPTIVE
    if(adaptive){
      for(it in 1:burnin){
        for(k in 1:reps){
          #local parameters
          aux <- sample_groups_vect_update_sigma_barker(theta, Y, X, tau, sigma,theta_mean, diag_sd, it)
          theta <- aux[[1]]
          sigma <- aux[[2]]
          theta_mean <- aux[[3]]
          diag_sd <- aux[[4]]
        }
        
        #sample tau
        tau <- sample_global_tau(theta, alpha1, alpha2)
      }
      if(verbose){
        print(sigma)
        print(c(min(diag_sd), max(diag_sd)))
      }
    }else{#not adaptive
      for(it in 1:burnin){
        for(k in 1:reps){
          #sample group parameters
          theta <- sample_groups_vect_barker(theta, Y, X, tau, sigma,theta_mean, diag_sd)
        }
        #sample tau
        tau <- sample_global_tau(theta, alpha1, alpha2)
      }
    }
    #print(paste("After burnin:",tau))
    #iterations
    Tau <- rep(0, iterations)
    Theta <- matrix(0, nrow = iterations, ncol = p)
    for(it in 1:iterations){
      for(k in 1:reps){
        #sample group parameters
        theta <- sample_groups_vect_barker(theta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau(theta, alpha1, alpha2)
      
      #save values
      Tau[it] <- tau
      Theta[it,] <- theta
    }
  }
  
  #RWM
  if(method == "rwm"){
    #warm start
    for(it in 1:(it_warm)){
      aux <- sample_groups_vect_update_sigma_rwm(theta, Y, X,tau, sigma,theta_mean, diag_sd, it)
      theta <- aux[[1]]
      sigma <- aux[[2]]
      sigma_store[i] <- sigma
      theta_mean <- aux[[3]]
      diag_sd <- aux[[4]]
    }
    if(verbose){
      print("Adaptively tuned parameters:")
      print(sigma)
      print(c(min(diag_sd), max(diag_sd)))
      plot(sigma_store, type = "l", main = "Adaptive sigma")
      cat("\n")
    }
    
    #burnin
    
    #IF ADAPTIVE
    if(adaptive){
      for(it in 1:burnin){
        for(k in 1:reps){
          #local parameters
          aux <- sample_groups_vect_update_sigma_barker(theta, Y, X, tau, sigma,theta_mean, diag_sd, it)
          theta <- aux[[1]]
          sigma <- aux[[2]]
          theta_mean <- aux[[3]]
          diag_sd <- aux[[4]]
        }
        
        #sample tau
        tau <- sample_global_tau(theta, alpha1, alpha2)
      }
      if(verbose){
        print(sigma)
        print(c(min(diag_sd), max(diag_sd)))
      }
    }
    else{
      for(it in 1:burnin){
        for(k in 1:reps){
          #sample group parameters
          theta <- sample_groups_vect_rwm(theta, Y, X, tau, sigma,theta_mean, diag_sd)
        }
        #sample tau
        tau <- sample_global_tau(theta, alpha1, alpha2)
      }
    }
    
    #iterations
    Tau <- rep(0, iterations)
    Theta <- matrix(0, nrow = iterations, ncol = p)
    for(it in 1:iterations){
      for(k in 1:reps){
        #sample group parameters
        theta <- sample_groups_vect_rwm(theta, Y, X, tau, sigma,theta_mean, diag_sd)
      }
      #sample tau
      tau <- sample_global_tau(theta, alpha1, alpha2)
      
      #save values
      Tau[it] <- tau
      Theta[it,] <- theta
    }
  }
  
  ESS_min <- c(ess(Tau),ess(rowSums2(Theta^2)), min(ess(Theta)), min(ess(Theta^2)))
  IAT_max <- iterations/ESS_min
  if(verbose){
    print("Iats for different parameters:")
    print(c(IAT_max))
    cat("\n")
  }
  
  to_return <- list()
  to_return[[1]] <- Tau
  to_return[[2]] <- Theta
  #to_return[[3]] <- sigma
  to_return[[4]] <- IAT_max
  
  return(to_return)
}

