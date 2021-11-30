############################################
# This is an R script used to fit the
# pooled Stan model and compute model
# evaluation diagnostics for the "Urinary
# biomarkers for pancreatic cancer" dataset.
############################################

library(rstan)
library(loo)

data_path <- "Kurssit/Bayesian Data Analysis/Project/Debernardi et al 2020 data.csv" #'~/your/path/to/dataset' # Replace with working path
data <- read.csv(data_path)

# Split the data into test subject groups
group_1 <- data[data$diagnosis == 1,]
group_2 <- data[data$diagnosis == 2,]
group_3 <- data[data$diagnosis == 3,]

# Organize input data for the Stan model 
stan_data <- list(
  N1 = length(group_1[,1]),
  N2 = length(group_2[,1]),
  N3 = length(group_3[,1]),
  y1 = group_1[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y2 = group_2[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y3 = group_3[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')]
)

# Fit the Stan model
pooled_fit <- stan(
  file = 'pooled.stan',
  data = stan_data
)

# Extract log-likelihoods for LOO evaluation
log_liks <- list(
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_1', merge_chains = FALSE),
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_2', merge_chains = FALSE),
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_3', merge_chains = FALSE)
)

# Initialize the diagnostics matrix
diagnostics <- matrix(0, nrow = 3, ncol = 2,
                      dimnames = list(
                        c('Group 1', 'Group2', 'Group 3'),
                        c('ELPD', 'P_eff')
                      )
)

# Compute diagnostics for each group
for (i in 1:3) {
  log_lik <- log_liks[[i]]
  r_eff <- relative_eff(exp(log_lik)) # Relative efficiency
  loo <- loo(log_lik, r_eff=r_eff) # LOO object
  
  estimates <- loo$estimates
  elpd <- estimates[1,1] # PSIS-LOO value
  p_eff <- estimates[2,1] # Effective number of parameters
  
  diagnostics[i,1] = elpd
  diagnostics[i,2] = p_eff
  print(loo, plot_k = TRUE)
}
