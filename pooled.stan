// -------------------------------------
// This is a pooled Stan model for the
// pancreatic cancer biomarker data set.  
// -------------------------------------

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
  vector<lower=0> [4] alpha;
  vector<lower=0> [4] beta;
}

model {
  // Priors
  for (j in 1:4) {
    alpha[j] ~ gamma(1, 1);
    beta[j] ~ gamma(0.5, 1);    
  }
  // Likelihoods
  for (j in 1:4) {
    y1[,j] ~ gamma(alpha[j], beta[j]);
    y2[,j] ~ gamma(alpha[j], beta[j]);
    y3[,j] ~ gamma(alpha[j], beta[j]);
  }
}

generated quantities {
  // Log-likelihoods of the posterior draws
  vector[4] log_lik_1[N1];
  vector[4] log_lik_2[N2];
  vector[4] log_lik_3[N3];
  // Group 1
  for (j in 1:4) {
    for (n in 1:N1) {
      log_lik_1[n,j] = gamma_lpdf(y1[n,j] | alpha[j], beta[j]);
    }
  }
  //Group 2
  for (j in 1:4) {
    for (n in 1:N2) {
      log_lik_2[n,j] = gamma_lpdf(y2[n,j] | alpha[j], beta[j]);
    }
  }
  //Group 3
  for (j in 1:4) {
    for (n in 1:N3) {
      log_lik_3[n,j] = gamma_lpdf(y3[n,j] | alpha[j], beta[j]);
    }
  }
}

