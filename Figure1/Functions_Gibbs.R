gibbs_lpost <- function(theta, y, m_data, mu, tau){
  res <- y*theta-m_data*log(1+exp(theta))+0.5*log(tau)-0.5*log(2*pi)
  res <- res-tau/2*(theta-mu)^2
  return(res)
}
gibbs_grad_lpost <- function(theta, y, m_data, mu, tau){
  res <- y-m_data*exp(theta)/(1+exp(theta))-tau*(theta-mu)
  return(res)
}
#SAMPLE GLOBAL PARAMETERS
sample_global_mu <- function(theta, tau){
  J <- length(theta)
  res <- mean(theta)+1/sqrt(J*tau)*rnorm(1)
  return(res)
}
sample_global_tau <- function(theta, mu, alpha, beta){
  J <- length(theta)
  res <- rgamma(1, alpha+J/2, beta+0.5*sum((theta-mean(theta))^2))
  return(res)
}
Gibbs_sampler <- function(Y, m, iterations, burnin, alpha, beta, mu0 = 0, tau0 = 1, verbose = F){
  #tau = precision of the groups
  J <- length(Y) #number of groups
  mu <- mu0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J)) #starting point
  tau <- tau0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J))
  
  #warm start
  theta <- mu+1/sqrt(tau)*rnorm(J) #group specific parameters
  it_warm <- as.integer(10*log(J))+1
  for(it in 1:it_warm){
    for(j in 1:J){
      theta[j] <- ars(1, gibbs_lpost, gibbs_grad_lpost, y = Y[j], m_data = m, mu = mu, tau = tau)
    }
  }
  
  for(it in 1:burnin){
    #sample group parameters
    for(j in 1:J){
      theta[j] <- ars(1, gibbs_lpost, gibbs_grad_lpost, y = Y[j], m_data = m, mu = mu, tau = tau)
    }
    #sample tau
    tau <- sample_global_tau(theta, mu, alpha, beta)
    #sample mu
    mu <- sample_global_mu(theta, tau)
  }
  
  Mu <- rep(0, iterations)
  Tau <- rep(0, iterations)
  Theta <- matrix(0, nrow = iterations, ncol = J)
  for(it in 1:iterations){
    #sample group parameters
    for(j in 1:J){
      theta[j] <- ars(1, gibbs_lpost, gibbs_grad_lpost, y = Y[j], m_data = m, mu = mu, tau = tau)
    }
    #sample tau
    tau <- sample_global_tau(theta, mu, alpha, beta)
    #sample mu
    mu <- sample_global_mu(theta, tau)
    
    #save values
    Tau[it] <- tau
    Mu[it] <- mu
    Theta[it,] <- theta
  }
  ESS_min <- min(c(ess(Mu), ess(Tau), ess(Theta)))
  IAT_max <- iterations/ESS_min
  
  to_return <- list()
  to_return[[1]] <- Mu
  to_return[[2]] <- Tau
  to_return[[3]] <- Theta
  to_return[[5]] <- IAT_max
  
  return(to_return)
}
