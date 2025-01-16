// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> m; //the number of observations per group
  int<lower=1> J; //the number of groups
  int<lower = 0, upper = m> y[J]; //the response variable
  real<lower=0> alpha; //first parameter gamma prior on tau
  real<lower=0> beta; //second parameter gamma prior on tau
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[J] theta; //local parameter
  real<lower=0> tau; //prior precision of theta
  real mu; //prior mean of theta
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  //priors global parameters
  tau ~ gamma(alpha, beta);
  mu ~ normal(0, sqrt(1000/tau));
  
  // priors local parameters
  for(j in 1:J){
    theta[j] ~ normal(mu, 1/sqrt(tau));
  }
  
  //likelihood
  for(j in 1:J){
    y[j] ~ binomial_logit(m, theta[j]);
  }
}

