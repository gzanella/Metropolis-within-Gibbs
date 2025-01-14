#computing logposterior
lpost_alias <- function(theta, parms){
  return(-lpost(parms[[1]],parms[[2]], theta, parms[[3]], parms[[4]]))
}
lpost <- function(y,x, theta, mu, tau){
  prod_xtheta <- x%*%theta
  
  res <- sum(y*prod_xtheta-log(1+exp(prod_xtheta))) #likelihood
  res <- res-sum(tau/2*(theta-mu)^2) #prior of theta
  res <- as.vector(res)
  
  return(res)
}
grad_lpost <- function(y,x,theta, mu, tau){
  prod_xtheta <- x%*%theta
  
  res <- matrix(y-exp(prod_xtheta)/(1+exp(prod_xtheta)), nrow = 1)%*%x #likelihood
  res <- res-tau*(theta-mu) #prior on theta
  res <- as.vector(res)
  
  return(res)
}

#GLOBAL PARAMETERS
sample_global_mu <- function(theta, tau){
  J <- nrow(theta)
  p <- ncol(theta)
  res <- colMeans2(theta)+1/sqrt(J*tau)*rnorm(p)
  
  return(res)
}
sample_global_tau <- function(theta, mu, alpha, beta){
  J <- nrow(theta)
  p <- ncol(theta)
  res <- rgamma(p, alpha+J/2, beta+0.5*(colSums2(theta^2)-J*colMeans2(theta)^2))
  
  return(res)
}

#BARKER
rbarker_diag<-function(theta,c,sigma,diag_sd){
  z<-sigma*rnorm(n=length(c),mean = 0,sd = diag_sd)
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  
  return(theta+z*b)
}
log_q_ratio_diag_barker<-function(x,y,grad_x, grad_y){
  beta1<-  c(-grad_y*(x-y))
  beta2<-  c(-grad_x*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
sample_groups_vect_barker <- function(theta, y, x, mu, tau, sigma, diag_sd){
  
  #proposal
  grad_theta <- grad_lpost(y,x,theta, mu, tau)
  
  barker <- rbarker_diag(theta,grad_theta,sigma, diag_sd)
  grad_barker <- grad_lpost(y,x,barker, mu, tau)
  
  #acceptance probability
  log_acc <- lpost(y,x,barker, mu, tau)-lpost(y,x,theta, mu, tau)+
    log_q_ratio_diag_barker(theta,barker,grad_theta,grad_barker)
  
  if(log(runif(1)) <= log_acc){
    theta <- barker
    #print("Accepted")
  }
  
  return(theta)
}
sample_groups_vect_update_sigma_barker <- function(theta, y,x, mu, tau, sigma, theta_mean, diag_sd, t){
  #proposal
  grad_theta <- grad_lpost(y,x,theta, mu, tau)
  
  barker <- rbarker_diag(theta,grad_theta,sigma, diag_sd)
  grad_barker <- grad_lpost(y,x,barker, mu, tau)
  
  
  #acceptance probability
  log_acc <- lpost(y,x,barker, mu, tau)-lpost(y,x,theta, mu, tau)+
    log_q_ratio_diag_barker(theta,barker,grad_theta,grad_barker)
  
  #print(theta_mean)
  #print(sigma)
  #print(diag_sd)
  #print(theta)
  #print(barker)
  #cat("\n")
  
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

#INDEPENDENT MH
sample_groups_vect_independent <- function(theta, y, x, mu, tau, reps){
  #compute modes
  parms <- list(y,x,mu,tau)
  means <- optim(par = mu, fn = lpost_alias,gr = NULL, parms, method = "BFGS")$par
  #means <- mu
  
  p <- length(theta)
  for(k in 1:reps){
    #sample proposal
    proposal <- means + 1/sqrt(tau)*rnorm(p)
    
    #compute alpha
    log_acc <- lpost(y,x,proposal, mu, tau)-lpost(y,x,theta, mu, tau)+
      sum(dnorm(theta, mean = means, sd = 1/sqrt(tau), log = T))-sum(dnorm(proposal, mean = means, sd = 1/sqrt(tau), log = T))
    
    #accept or reject
    if(log(runif(1)) <= log_acc){
      theta <- proposal
      #print("Accepted")
    }
  }
  
  return(theta)
}

#Random Walk MH
sample_groups_vect_rwm <- function(theta, y, x, mu, tau, sigma, diag_sd){
  
  p <- length(theta)
  #proposal
  proposal <- theta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost(y,x,proposal, mu, tau)-lpost(y,x,theta, mu, tau)
  
  if(log(runif(1)) <= log_acc){
    theta <- proposal
    #print("Accepted")
  }
  
  return(theta)
}
sample_groups_vect_update_sigma_rwm <- function(theta, y,x, mu, tau, sigma, theta_mean, diag_sd, t){
  #proposal
  p <- length(theta)
  proposal <- theta+sigma*rnorm(p, 0, diag_sd)
  
  #acceptance probability
  log_acc <- lpost(y,x,proposal, mu, tau)-lpost(y,x,theta, mu, tau)
  
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
MH_sampler_vect <- function(Y,X, m, iterations, burnin, alpha, beta, sigma = 2.4, mu0 = 0, tau0 = 1, adaptive = F, reps = 1, verbose = F, method = "barker"){
  #Y = list: each element is a m-vector
  #X = list: each element is a matrix
  #m = number of observations per group
  #alpha, beta = prior parameters for tau
  #sigma = stepsize barker/rwm
  #mu0, tau0 = true parameters
  #reps = number of iterations MH per group
  #method = barker/independent/rwm
  
  J <- length(Y) #number of groups
  p <- ncol(X[[1]]) #number of variables
  
  #starting point
  mu <- mu0 + runif(p, min = -1/sqrt(J), max = 1/sqrt(J))
  tau <- tau0 + runif(p, min = -1/sqrt(J), max = 1/sqrt(J))
  it_warm <- as.integer(100*log(J))+1 #number iterations warm start
  
  #preconditioning
  sigma <- sigma
  theta_mean <- rep(0, p)
  diag_sd <- rep(1, p)
  
  #BARKER
  theta <- matrix(0, nrow = J, ncol = p)
  if(method == "barker"){
    #warm start
    for(j in 1:J){
      theta[j,] <- mu+1/sqrt(tau)*rnorm(p) #group specific parameters
      for(it in 1:(it_warm*reps)){
        aux <- sample_groups_vect_update_sigma_barker(theta[j,], Y[[j]], X[[j]], mu, tau, sigma,theta_mean, diag_sd, it)
        theta[j,] <- aux[[1]]
        sigma <- aux[[2]]
        theta_mean <- aux[[3]]
        diag_sd <- aux[[4]]
      }
    }
    
    #burnin
    
    #IF ADAPTIVE
    if(adaptive){
      for(it in 1:(reps*burnin)){
        #sample group parameters
        for(j in 1:J){
          aux <- sample_groups_vect_update_sigma_barker(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, theta_mean, diag_sd, it+it_warm)
          theta[j,] <- aux[[1]]
          sigma <- aux[[2]]
          theta_mean <- aux[[3]]
          diag_sd <- aux[[4]]
        }
        #sample tau
        tau <- sample_global_tau(theta, mu, alpha, beta)
        #sample mu
        mu <- sample_global_mu(theta, tau)
      }
      if(verbose){
        print(sigma)
        print(diag_sd)
      }
    }
    else{
      for(it in 1:burnin){
        for(j in 1:J){
          for(k in 1:reps){
            #sample group parameters
            theta[j,] <- sample_groups_vect_barker(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, diag_sd)
          }
        }
        #sample tau
        tau <- sample_global_tau(theta, mu, alpha, beta)
        #sample mu
        mu <- sample_global_mu(theta, tau)
      }
    }
    
    #iterations
    Mu <- matrix(0, nrow = iterations, ncol = p)
    Tau <- matrix(0, nrow = iterations, ncol = p)
    Theta <- vector(mode='list', length=p)
    for(k in 1:p){
      Theta[[k]] <- matrix(0, nrow = iterations, ncol = J)
    }
    for(it in 1:iterations){
      #sample group parameters
      for(k in 1:reps){
        for(j in 1:J){
          theta[j,] <- sample_groups_vect_barker(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, diag_sd)
        }
      }
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
      
      #save values
      Tau[it,] <- tau
      Mu[it,] <- mu
      for(k in 1:p){
        Theta[[k]][it,] <- theta[,k]
      }
    }
    ESS_min <- min(c(ess(Mu), ess(Tau), sapply(Theta, ess)))
    IAT_max <- iterations/ESS_min
    if(verbose){
      print(c(iterations/c(ess(Mu), ess(Tau)), IAT_max))
    }
  }
  
  #INDEPENDENT
  theta <- matrix(0, nrow = J, ncol = p)
  if(method == "independent"){
    #warm start
    for(j in 1:J){
      theta[j,] <- mu+1/sqrt(tau)*rnorm(p) #group specific parameters
      for(it in 1:it_warm){
        theta[j,] <- sample_groups_vect_independent(theta[j,], Y[[j]], X[[j]], mu, tau, reps)
      }
    }
    
    #burnin
    
  
    for(it in 1:burnin){
      for(j in 1:J){
        #sample group parameters
        theta[j,] <- sample_groups_vect_independent(theta[j,], Y[[j]], X[[j]], mu, tau, reps)
      }
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
      }
    
    #iterations
    Mu <- matrix(0, nrow = iterations, ncol = p)
    Tau <- matrix(0, nrow = iterations, ncol = p)
    Theta <- vector(mode='list', length=p)
    for(k in 1:p){
      Theta[[k]] <- matrix(0, nrow = iterations, ncol = J)
    }
    for(it in 1:iterations){
      #sample group parameters
      for(j in 1:J){
        theta[j,] <- sample_groups_vect_independent(theta[j,], Y[[j]], X[[j]], mu, tau, reps)
      }
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
      
      #save values
      Tau[it,] <- tau
      Mu[it,] <- mu
      for(k in 1:p){
        Theta[[k]][it,] <- theta[,k]
      }
    }
    ESS_min <- min(c(ess(Mu), ess(Tau), sapply(Theta, ess)))
    IAT_max <- iterations/ESS_min
    if(verbose){
      print(c(iterations/c(ess(Mu), ess(Tau)), IAT_max))
    }
  }
  
  #RWM
  theta <- matrix(0, nrow = J, ncol = p)
  if(method == "rwm"){
    #warm start
    for(j in 1:J){
      theta[j,] <- mu+1/sqrt(tau)*rnorm(p) #group specific parameters
      for(it in 1:(it_warm*reps)){
        aux <- sample_groups_vect_update_sigma_rwm(theta[j,], Y[[j]], X[[j]], mu, tau, sigma,theta_mean, diag_sd, it)
        theta[j,] <- aux[[1]]
        sigma <- aux[[2]]
        theta_mean <- aux[[3]]
        diag_sd <- aux[[4]]
      }
    }
    
    #burnin
    
    #IF ADAPTIVE
    if(adaptive){
      for(it in 1:(reps*burnin)){
        #sample group parameters
        for(j in 1:J){
          aux <- sample_groups_vect_update_sigma_rwm(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, theta_mean, diag_sd, it+it_warm)
          theta[j,] <- aux[[1]]
          sigma <- aux[[2]]
          theta_mean <- aux[[3]]
          diag_sd <- aux[[4]]
        }
        #sample tau
        tau <- sample_global_tau(theta, mu, alpha, beta)
        #sample mu
        mu <- sample_global_mu(theta, tau)
      }
      if(verbose){
        print(sigma)
        print(diag_sd)
      }
    }
    else{
      for(it in 1:burnin){
        for(j in 1:J){
          for(k in 1:reps){
            #sample group parameters
            theta[j,] <- sample_groups_vect_rwm(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, diag_sd)
          }
        }
        #sample tau
        tau <- sample_global_tau(theta, mu, alpha, beta)
        #sample mu
        mu <- sample_global_mu(theta, tau)
      }
    }
    
    #iterations
    Mu <- matrix(0, nrow = iterations, ncol = p)
    Tau <- matrix(0, nrow = iterations, ncol = p)
    Theta <- vector(mode='list', length=p)
    for(k in 1:p){
      Theta[[k]] <- matrix(0, nrow = iterations, ncol = J)
    }
    for(it in 1:iterations){
      #sample group parameters
      for(k in 1:reps){
        for(j in 1:J){
          theta[j,] <- sample_groups_vect_rwm(theta[j,], Y[[j]], X[[j]], mu, tau, sigma, diag_sd)
        }
      }
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
      
      #save values
      Tau[it,] <- tau
      Mu[it,] <- mu
      for(k in 1:p){
        Theta[[k]][it,] <- theta[,k]
      }
    }
    ESS_min <- min(c(ess(Mu), ess(Tau), sapply(Theta, ess)))
    IAT_max <- iterations/ESS_min
    if(verbose){
      print(c(iterations/c(ess(Mu), ess(Tau)), IAT_max))
    }
  }
  
  to_return <- list()
  to_return[[1]] <- Mu
  to_return[[2]] <- Tau
  to_return[[3]] <- Theta
  to_return[[4]] <- sigma
  to_return[[5]] <- IAT_max
  
  return(to_return)
}