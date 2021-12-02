// ------------------------------------------
// This is a pooled Stan model for the
// "Urinary biomarkers for pancreatic cancer"
// data set.  
// ------------------------------------------

data {
  // Lengths of the observation vectors
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N3;
  // Observations of the 4 protein levels for each group
  vector[4] y1[N1];
  vector[4] y2[N2];
  vector[4] y3[N3];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector [3] mu;
  vector<lower=0> [3] sigma;
}

model {
  // Priors
  alpha ~ gamma(1, 1);
  beta ~ gamma(0.5, 1);  
  for (j in 1:3) {
    mu[j] ~ normal(0,20);
    sigma[j] ~ gamma(1, 1);
  }
  // Likelihoods
  y1[,1] ~ gamma(alpha, beta);
  y2[,1] ~ gamma(alpha, beta);
  y3[,1] ~ gamma(alpha, beta);
  for (j in 1:3) {
    log(y1[,j+1]) ~ normal(mu[j], sigma[j]);
    log(y2[,j+1]) ~ normal(mu[j], sigma[j]);
    log(y3[,j+1]) ~ normal(mu[j], sigma[j]);
  }
}

generated quantities {
  // Posterior predictive distributions
  vector[4] ypred_1;
  vector[4] ypred_2;
  vector[4] ypred_3;
  // Log-likelihoods of the posterior draws
  vector[4] log_lik_1[N1];
  vector[4] log_lik_2[N2];
  vector[4] log_lik_3[N3];
  
  // Group 1
  ypred_1[1] = gamma_rng(alpha, beta);
  for (n in 1:N1)
    log_lik_1[n,1] = gamma_lpdf(y1[n,1] | alpha, beta);
    
  for (j in 1:3) {
    ypred_1[j+1] = normal_rng(mu[j], sigma[j]);
    for (n in 1:N1) {
      log_lik_1[n,j+1] = normal_lpdf(log(y1[n,j+1]) | mu[j], sigma[j]);
    }
  }
  //Group 2
  ypred_2[1] = gamma_rng(alpha, beta);
  for (n in 1:N2)
    log_lik_2[n,1] = gamma_lpdf(y2[n,1] | alpha, beta);
    
  for (j in 1:3) {
    ypred_2[j+1] = normal_rng(mu[j], sigma[j]);
    for (n in 1:N2) {
      log_lik_2[n,j+1] = normal_lpdf(log(y2[n,j+1]) | mu[j], sigma[j]);
    }
  }
  //Group 3
  ypred_3[1] = gamma_rng(alpha, beta);
  for (n in 1:N3)
      log_lik_3[n,1] = gamma_lpdf(y3[n,1] | alpha, beta);
      
  for (j in 1:3) {
    ypred_3[j+1] = normal_rng(mu[j], sigma[j]);
    for (n in 1:N3) {
      log_lik_3[n,j+1] = normal_lpdf(log(y3[n,j+1]) | mu[j], sigma[j]);
    }
  }
}

