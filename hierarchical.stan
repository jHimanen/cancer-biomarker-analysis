// This is a hierarchical Stan model for the
// pancreatic cancer biomarker data set.  

// The input data
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

transformed data {
  vector[3] log_y1[N1];
  vector[3] log_y2[N2];
  vector[3] log_y3[N3];
  
  for(j in 1:3){
    log_y1[,j] = log(y1[,j+1]);
    log_y2[,j] = log(y2[,j+1]);
    log_y3[,j] = log(y3[,j+1]);
  }
}

// The parameters accepted by the model
parameters {
  // Hyperparamters
  vector<lower=0>[4] shape[2];
  vector<lower=0>[4] scale[2];
  // Parameters
  vector<lower=0>[4] alpha[3];
  vector<lower=0>[4] beta[3];
}

// The model to be estimated
model {
  //k->number of parameters in priors
  //j->number of proteins
  // hyperpriors
  for (k in 1:2){
    shape[k,] ~ gamma(1,1);
    scale[k,] ~ gamma(1,1);
  }
  //priors
  for (j in 1:4){
    alpha[,j] ~ normal(shape[1,j], scale[1,j]);
    beta[,j] ~ gamma(shape[2,j], scale[2,j]);
  }
  //Likelihood
  y1[,1] ~ gamma(alpha[1,1], beta[1,1]);
  y2[,1] ~ gamma(alpha[2,1], beta[2,1]);
  y3[,1] ~ gamma(alpha[3,1], beta[3,1]);
  for (j in 1:3){
    log_y1[,j] ~ normal(alpha[1,j+1], beta[1,j+1]);
    log_y2[,j] ~ normal(alpha[2,j+1], beta[2,j+1]);
    log_y3[,j] ~ normal(alpha[3,j+1], beta[3,j+1]);
  }
}

