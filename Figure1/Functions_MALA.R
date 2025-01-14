logpost <- function(y,m,theta, mu, tau, alpha, beta){
  J <- length(theta)
  #likelihood
  res <- sum(y*theta)-m*sum(log(1+exp(theta)))
  #local parameters
  res <- res + J/2*log(tau)-tau/2*sum((theta-mu)^2)
  #global parameters
  res <- res + (alpha-1)*log(tau)-beta*tau
  
  return(res)
}
grad_theta <- function(y, m, theta, mu, tau){
  res <- y-m*exp(theta)/(1+exp(theta))-tau*(theta-mu)
  return(res)
}
grad_mu <- function(theta, mu, tau){
  res <- tau*sum(theta-mu)
  return(res)
}
grad_tau <- function(theta, mu, tau, alpha, beta){
  J <- length(theta)
  res <- (J/2+alpha-1)/tau-0.5*sum((theta-mu)^2)-beta
  return(res)
}
MALA_sampler <- function(Y, m, iterations, burnin, alpha, beta, gamma0,sd, mu0 = 0, tau0 = 1, reps, verbose = F){
  #tau = precision of the groups
  J <- length(Y) #number of groups
  mu <- mu0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J)) #starting point
  tau <- tau0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J))
  logtau <- log(tau)
  
  #warm start
  theta <- mu+1/sqrt(tau)*rnorm(J) #group specific parameters
  it_warm <- as.integer(10*log(J))+1
  for(it in 1:it_warm){
    theta <- sample_groups_vect(theta, Y, mu, tau, 1)
  }
  #burnin
  Gamma0 <- rep(0, reps*burnin)
  for(it in 1:(reps*burnin)){
    #stepsize
    gamma <- gamma0*sd/(J+2)^(1/6)
    #values old
    log_post <- logpost(Y,m,theta, mu, tau, alpha, beta)
    grad <- c(grad_mu(theta, mu, tau), 
              grad_tau(theta, mu, tau, alpha, beta),
              grad_theta(Y, m, theta, mu, tau))
    #proposal
    mu_new <- mu + gamma[1]^2*grad[1]/2+gamma[1]*rnorm(1)
    logtau_new <- logtau + gamma[2]^2*grad[2]/2+gamma[2]*rnorm(1)
    tau_new <- exp(logtau_new)
    theta_new <- theta + gamma[3:(J+2)]^2*grad[3:(J+2)]/2+gamma[3:(J+2)]*rnorm(J)
    #values new
    log_post_new <- logpost(Y,m,theta_new, mu_new, tau_new, alpha, beta)
    grad_new <- c(grad_mu(theta_new, mu_new, tau_new), 
              grad_tau(theta_new, mu_new, tau_new, alpha, beta),
              grad_theta(Y, m, theta_new, mu_new, tau_new))
    #log acceptance
      #likelihood
    log_acc <- log_post_new-log_post
      #mu
    log_acc <- log_acc+dnorm(mu, mean = mu_new+gamma[1]^2*grad_new[1]/2, sd = gamma[1], log = T)-dnorm(mu_new, mean = mu+gamma[1]^2*grad[1]/2, sd = gamma[1], log = T)
      #tau
    log_acc <- log_acc+dnorm(logtau, mean = logtau_new+gamma[2]^2*grad_new[2]/2, sd = gamma[2], log = T)-dnorm(logtau_new, mean = logtau+gamma[2]^2*grad[2]/2, sd = gamma[2], log = T)
      #theta
    log_acc <- log_acc+sum(dnorm(theta, mean = theta_new+gamma[3:(J+2)]^2*grad_new[3:(J+2)]/2, sd = gamma[3:(J+2)], log = T)-dnorm(theta_new, mean = theta+gamma[3:(J+2)]^2*grad[3:(J+2)]/2, sd = gamma[3:(J+2)], log = T))
    if(log(runif(1)) <= log_acc){
      theta <- theta_new
      mu <- mu_new
      tau <- tau_new
      logtau <- logtau_new
    }
    #update gamma0
    acc <- pmin(exp(log_acc), 1)
    gamma0_log <- log(gamma0)+it^(-1)*(acc-0.57)
    gamma0 <- exp(gamma0_log)
    Gamma0[it] <- gamma0
  }
  if(verbose){
    plot(Gamma0, type = "l")
  }
  Mu <- rep(0, iterations)
  Tau <- rep(0, iterations)
  Theta <- matrix(0, nrow = iterations, ncol = J)
  n_acc <- 0
  for(it in 1:iterations){
    for(j in 1:reps){
      #values old
      log_post <- logpost(Y,m,theta, mu, tau, alpha, beta)
      grad <- c(grad_mu(theta, mu, tau), 
                grad_tau(theta, mu, tau, alpha, beta),
                grad_theta(Y, m, theta, mu, tau))
      #proposal
      mu_new <- mu + gamma[1]^2*grad[1]/2+gamma[1]*rnorm(1)
      logtau_new <- logtau + gamma[2]^2*grad[2]/2+gamma[2]*rnorm(1)
      tau_new <- exp(logtau_new)
      theta_new <- theta + gamma[3:(J+2)]^2*grad[3:(J+2)]/2+gamma[3:(J+2)]*rnorm(J)
      #values new
      log_post_new <- logpost(Y,m,theta_new, mu_new, tau_new, alpha, beta)
      grad_new <- c(grad_mu(theta_new, mu_new, tau_new), 
                    grad_tau(theta_new, mu_new, tau_new, alpha, beta),
                    grad_theta(Y, m, theta_new, mu_new, tau_new))
      #log acceptance
      #likelihood
      log_acc <- log_post_new-log_post
      #mu
      log_acc <- log_acc+dnorm(mu, mean = mu_new+gamma[1]^2*grad_new[1]/2, sd = gamma[1], log = T)-dnorm(mu_new, mean = mu+gamma[1]^2*grad[1]/2, sd = gamma[1], log = T)
      #tau
      log_acc <- log_acc+dnorm(logtau, mean = logtau_new+gamma[2]^2*grad_new[2]/2, sd = gamma[2], log = T)-dnorm(logtau_new, mean = logtau+gamma[2]^2*grad[2]/2, sd = gamma[2], log = T)
      #theta
      log_acc <- log_acc+sum(dnorm(theta, mean = theta_new+gamma[3:(J+2)]^2*grad_new[3:(J+2)]/2, sd = gamma[3:(J+2)], log = T)-dnorm(theta_new, mean = theta+gamma[3:(J+2)]^2*grad[3:(J+2)]/2, sd = gamma[3:(J+2)], log = T))
      if(log(runif(1)) <= log_acc){
        n_acc <- n_acc+1
        theta <- theta_new
        mu <- mu_new
        tau <- tau_new
        logtau <- logtau_new
      }
    }
    #save values
    Tau[it] <- tau
    Mu[it] <- mu
    Theta[it,] <- theta
  }
  #print(c(ess(Mu), ess(Tau),ess(Theta)))
  ESS_min <- min(c(ess(Mu), ess(Tau), ess(Theta)))
  #print(ESS_min)
  IAT_max <- iterations/ESS_min
  if(verbose){
    print(paste("Proportion of accepted values:",n_acc/(reps*iterations)))
    print(summary(iterations/ess(Theta)))
    print(c(iterations/c(ess(Mu), ess(Tau)), IAT_max))
  }
  
  return(list(mu= Mu, tau = Tau, IAT_max = IAT_max))
}
