// ------------------------------------------
// This is a hierarchical Stan model for the
// "Urinary biomarkers for pancreatic cancer"
// dataset.  
// ------------------------------------------

data {
  // Lengths of the observation vectors
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N3;
  // Observations for each group
  vector[4] y1[N1];
  vector[4] y2[N2];
  vector[4] y3[N3];
}

parameters {
  // Hyperparamters
  real<lower=0> alphaP[2];
  real<lower=0> betaP[2];
  vector[3] muP[2];
  vector<lower=0>[3] sigmaP[2];
  // Parameters
  real<lower=0> alpha[3];
  real<lower=0> beta[3];
  vector [3] mu[3];
  vector<lower=0> [3] sigma[3];
}

model {
  // k-> number of parameters in priors
  // j-> number of proteins
  // Hyperpriors
  for (k in 1:2){
    alphaP[k] ~ gamma(1,1);
    betaP[k] ~ gamma(1,1);
    sigmaP[k,] ~ gamma(1,1);
  }
  muP[1,] ~ normal(0,20);
  muP[2,] ~ gamma(1,1);
  // Priors
  alpha ~ gamma(alphaP[1], betaP[1]);
  beta ~ gamma(alphaP[2], betaP[2]);
  for (j in 1:3){
    mu[j] ~ normal(muP[1,j], sigmaP[1,j]);
    sigma[j] ~ gamma(muP[2,j], sigmaP[2,j]);
  }
  // Likelihood
  y1[,1] ~ gamma(alpha[1], beta[1]);
  y2[,1] ~ gamma(alpha[2], beta[2]);
  y3[,1] ~ gamma(alpha[3], beta[3]);
  for (j in 1:3){
    log(y1[,j+1]) ~ normal(mu[1,j], sigma[1,j]);
    log(y2[,j+1]) ~ normal(mu[2,j], sigma[2,j]);
    log(y3[,j+1]) ~ normal(mu[3,j], sigma[3,j]);
  }
}

