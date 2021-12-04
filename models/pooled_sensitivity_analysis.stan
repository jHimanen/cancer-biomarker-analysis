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
  vector[8] p;
}

transformed data {
  vector[3] log_y1[N1];
  vector[3] log_y2[N2];
  vector[3] log_y3[N3];
  
  for (j in 1:3){
    log_y1[,j] = log(y1[,j+1]);
    log_y2[,j] = log(y2[,j+1]);
    log_y3[,j] = log(y3[,j+1]);
  }
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector [3] mu;
  vector<lower=0> [3] sigma;
}

model {
  // Priors
  alpha ~ gamma(p[1], p[2]);
  beta ~ gamma(p[3], p[4]);  
  for (j in 1:3) {
    mu[j] ~ normal(p[5], p[6]);
    sigma[j] ~ gamma(p[7], p[8]);
  }
  // Likelihoods
  y1[,1] ~ gamma(alpha, beta);
  y2[,1] ~ gamma(alpha, beta);
  y3[,1] ~ gamma(alpha, beta);
  for (j in 1:3) {
    log_y1[,j] ~ normal(mu[j], sigma[j]);
    log_y2[,j] ~ normal(mu[j], sigma[j]);
    log_y3[,j] ~ normal(mu[j], sigma[j]);
  }
}
