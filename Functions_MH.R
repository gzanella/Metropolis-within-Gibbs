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
#BARKER VECTORIZED
log_q_ratio_barker_vect<-function(x,y,grad_x,grad_y){
  # x: current location (vector)
  # y: proposed location (vector)
  # grad_x: target log-posterior gradient at x (vector)
  # grad_y: target log-posterior gradient at y (vector)
  beta1<-  c(-grad_y*(x-y))
  beta2<-  c(-grad_x*(y-x))
  return(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2)))
  ))
}
#SAMPLE LOCAL PARAMETERS VECTORIZED
lpost <- function(y,m,theta, mu, tau){
  res <- y*theta-m*log(1+exp(theta))+0.5*log(tau)-0.5*log(2*pi)
  res <- res-tau/2*(theta-mu)^2
  return(res)
}
grad_lpost <- function(y,m,theta, mu, tau){
  res <- y-m*exp(theta)/(1+exp(theta))-tau*(theta-mu)
  return(res)
}
rbarker <-function(x,grad,sigma){
  # x: current location (vector)
  # grad: target log-posterior gradient (vector)
  # sigma: proposal stepsize (scalar)
  z<-rnorm(n=length(grad),mean = sigma,sd = 0.1*sigma)
  b<-2*(runif(n=length(grad))< 1/(1+exp(-grad*z)))-1
  return(x+z*b)
}
log_q_ratio_barker<-function(x,y,grad_x,grad_y){
  # x: current location (vector)
  # y: proposed location (vector)
  # grad_x: target log-posterior gradient at x (vector)
  # grad_y: target log-posterior gradient at y (vector)
  beta1<-  c(-grad_y*(x-y))
  beta2<-  c(-grad_x*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
sample_groups_vect <- function(theta, y, mu, tau, sigma){
  #non adaptive case
  J <- length(theta)
  grad_theta <- grad_lpost(y,m,theta, mu, tau)
  barker <- rbarker(theta,grad_theta,sigma)
  grad_barker <- grad_lpost(y,m,barker, mu, tau)
  log_acc <- lpost(y,m,barker, mu, tau)-lpost(y,m,theta, mu, tau)+
    log_q_ratio_barker_vect(theta,barker,grad_theta,grad_barker)
  U <- runif(J)
  theta[log(U) <= log_acc] <- barker[log(U) <= log_acc]
  return(theta)
}
sample_groups_vect_update_sigma <- function(theta, y, mu, tau, sigma, t){
  #adaptive case
  J <- length(theta)
  grad_theta <- grad_lpost(y,m,theta, mu, tau)
  barker <- rbarker(theta,grad_theta,sigma)
  grad_barker <- grad_lpost(y,m,barker, mu, tau)
  log_acc <- lpost(y,m,barker, mu, tau)-lpost(y,m,theta, mu, tau)+
    log_q_ratio_barker_vect(theta,barker,grad_theta,grad_barker)
  U <- runif(J)
  theta[log(U) <= log_acc] <- barker[log(U) <= log_acc]
  #update sigma
  acc <- pmin(1, exp(log_acc))
  sigma_log <- log(sigma)+t^(-0.8)*(mean(acc)-0.4)
  sigma <- exp(sigma_log)
  #print(summary(sigma))
  to_return <- list(theta, sigma)#matrix(c(theta, sigma), ncol = 2, byrow = F)
  return(to_return)
}
#METROPOLIS HASTINGS VECTORIZED
MH_sampler_vect <- function(Y, m, iterations, burnin, alpha, beta, sigma = 2.4, mu0 = 0, tau0 = 1, adaptive = F, reps = 1, verbose = F){
  #tau = precision of the groups
  J <- length(Y) #number of groups
  #if(adaptive){
   # sigma <- rep(sigma, J) #sigma Barker proposal
  #}
  mu <- mu0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J)) #starting point
  tau <- tau0 + runif(1, min = -1/sqrt(J), max = 1/sqrt(J))
  
  #warm start
  theta <- mu+1/sqrt(tau)*rnorm(J) #group specific parameters
  it_warm <- as.integer(10*log(J))+1
  for(it in 1:it_warm){
    theta <- sample_groups_vect(theta, Y, mu, tau, sigma)
  }
  #IF ADAPTIVE
  Sigma <- rep(0, burnin)
  if(adaptive){
    for(it in 1:(reps*burnin)){
      #sample group parameters
      aux <- sample_groups_vect_update_sigma(theta, Y, mu, tau, sigma, it)
      theta <- aux[[1]]
      sigma <- aux[[2]]
      Sigma[it] <- sigma
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
    }
    plot(Sigma, type = "l", main = "Sigma")
  }
  else{
    for(it in 1:burnin){
      for(j in 1:reps){
        #sample group parameters
        theta <- sample_groups_vect(theta, Y, mu, tau, sigma)
      }
      #sample tau
      tau <- sample_global_tau(theta, mu, alpha, beta)
      #sample mu
      mu <- sample_global_mu(theta, tau)
    }
  }
  Mu <- rep(0, iterations)
  Tau <- rep(0, iterations)
  Theta <- matrix(0, nrow = iterations, ncol = J)
  for(it in 1:iterations){
    for(j in 1:reps){
      #sample group parameters
      theta <- sample_groups_vect(theta, Y, mu, tau, sigma)
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
  if(verbose){
    print(c(iterations/c(ess(Mu), ess(Tau)), IAT_max))
  }
  
  to_return <- list()
  to_return[[1]] <- Mu
  to_return[[2]] <- Tau
  to_return[[3]] <- Theta
  to_return[[4]] <- sigma
  to_return[[5]] <- IAT_max
  
  return(to_return)
}


#HMC leapfrogs
leap_ESS<-function(fit,alg_stan,t_burn,print=FALSE){
  # compute num leap
  sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
  sampler_params_chain1 <- sampler_params[[1]]
  if(alg_stan=="NUTS"){
    leap_num<-sampler_params_chain1[,4]
  }
  if(alg_stan=="HMC"){
    leap_num<-sampler_params_chain1[,3]/sampler_params_chain1[,2]
  }
  if(print){
    print(paste("Mean number of leapfrog steps per iteration=",round(mean(leap_num))))
    print(paste("Mean number of leapfrog steps per iteration (warm-up)=",round(mean(leap_num[1:(length(leap_num)/2)]))))
  }
  # COMPUTE ESS WITH SAME PROCEDURE AS GIBBS
  samples<-extract(object=fit, permuted = FALSE, inc_warmup = TRUE,  include = TRUE)
  samples.stan<-samples[-c(1:t_burn),1,]
  ess_stan<-apply(samples.stan,2,ess)
  if(print){
    print("ESS:")
    print(summary(ess_stan))
    print(" ESS per grad ev:")
    print(summary(ess_stan/sum(leap_num)))
  }
  return(list(leap_num=mean(leap_num),ess_stan=ess_stan))
}