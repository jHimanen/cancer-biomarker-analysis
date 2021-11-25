// This is a pooled Stan model for the
// pancreatic cancer biomarker data set.  

// The input data
data {
  // Lengths of the observation vectors
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N3;
  // Observations for each group
  vector[N1] y1;
  vector[N2] y2;
  vector[N3] y3;
}

// The parameters accepted by the model
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

// The model to be estimated
model {
  //priors
  alpha ~ gamma(1, 1);
  beta ~ gamma(0.5, 1);
  //Likelihood
  y1 ~ gamma(alpha, beta);
  y2 ~ gamma(alpha, beta);
  y3 ~ gamma(alpha, beta);
}

