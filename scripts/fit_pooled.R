##############################################
# This is an R script used to fit the
# pooled Stan model and compute convergence, 
# model evaluation, and posterior predictive
# diagnostics for the "Urinary
# biomarkers for pancreatic cancer" dataset.
##############################################

library(rstan)
library(loo)
library(ggplot2)

data_path <- '.../Debernardi et al 2020 data.csv' # REPLACE WITH WORKING PATH
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
  file = '.../cancer-biomarker-analysis/models/pooled.stan', # REPLACE WITH WORKING PATH
  data = stan_data
)

# Store convergence diagnostics
results <- monitor(pooled_fit)
selected_res <- results[c(1:20), c('mean', 'sd', 'n_eff', 'Rhat', 'Q5', 'Q50', 'Q95')]
res_df <- as.data.frame(selected_res)
write.csv(res_df,'.../cancer-biomarker-analysis/diagnostic_data/pooled_res.csv') # REPLACE WITH WORKING PATH

# Extract log-likelihoods for LOO evaluation
log_liks <- list(
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_1', merge_chains = FALSE),
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_2', merge_chains = FALSE),
  extract_log_lik(pooled_fit, parameter_name = 'log_lik_3', merge_chains = FALSE)
)

# Initialize the diagnostics matrix
diagnostics <- matrix(0, nrow = 3, ncol = 2,
                      dimnames = list(
                        c('Group 1', 'Group 2', 'Group 3'),
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

pool_eval <- as.data.frame(diagnostics)
write.csv(pool_eval, '.../cancer-biomarker-analysis/diagnostic_data/pooled_eval.csv') # REPLACE WITH WORKING PATH

# Visualize posterior distributions

draws <- as.data.frame(pooled_fit)

pCreatinine <- ggplot() + ggtitle('Posterior distribution of creatinine')
x <- seq(0,5,0.05)
for (i in seq(1,4000,50)) {
  df <- data.frame(x=x,y=dgamma(x,draws$alpha[i],draws$beta[i]))
  pCreatinine <- pCreatinine + geom_line(data=df,aes(x,y),col='lightblue',alpha=0.4)
}
pCreatinine

pLYVE1 <- ggplot() + ggtitle('Posterior distribution of LYVE1')
x <- seq(-10,10,0.1)
for (i in seq(1,4000,50)) {
  df <- data.frame(x=x,y=dnorm(x,draws$`mu[1]`[i],draws$`sigma[1]`[i]))
  pLYVE1 <- pLYVE1 + geom_line(data=df,aes(x,y),col='lightblue',alpha=0.4)
}
pLYVE1

pREG1B <- ggplot() + ggtitle('Posterior distribution of REG1B')
x <- seq(-5,12,0.1)
for (i in seq(1,4000,50)) {
  df <- data.frame(x=x,y=dnorm(x,draws$`mu[2]`[i],draws$`sigma[2]`[i]))
  pREG1B <- pREG1B + geom_line(data=df,aes(x,y),col='lightblue',alpha=0.4)
}
pREG1B

pTFF1 <- ggplot() + ggtitle('Posterior distribution of TFF1')
x <- seq(-7,15,0.1)
for (i in seq(1,4000,50)) {
  df <- data.frame(x=x,y=dnorm(x,draws$`mu[3]`[i],draws$`sigma[3]`[i]))
  pTFF1 <- pTFF1 + geom_line(data=df,aes(x,y),col='lightblue',alpha=0.4)
}
pTFF1

# Posterior predictive check

# Visualizing posterior predictive distributions
# Note that ypred_i are from same distributions
stan_hist(pooled_fit, pars = c('ypred_1[1]','ypred_1[2]','ypred_1[3]','ypred_1[4]'))
stan_hist(pooled_fit, pars = c('ypred_2[1]','ypred_2[2]','ypred_2[3]','ypred_2[4]'))
stan_hist(pooled_fit, pars = c('ypred_3[1]','ypred_3[2]','ypred_3[3]','ypred_3[4]'))

# Drawing N1 samples from ypred1_[1] four times and plotting the histogram

pPostCheck_C <- ggplot() + geom_histogram(data=data,aes(x=creatinine),fill='white',color="black", binwidth = 0.1) +
  ggtitle("Replicated datasets of creatinine compared to the original data")
for (j in 1:4) {
  df <- data.frame(values=sample(draws$`ypred_1[1]`,(stan_data$N1+stan_data$N2+stan_data$N3)))
  pPostCheck_C <- pPostCheck_C + geom_histogram(data=df,aes(x=values),fill=j+1, alpha = 0.2, binwidth = 0.1)
}
pPostCheck_C 

pPostCheck_L <- ggplot() + geom_histogram(data=data,aes(x=log(LYVE1)),fill='white',color="black", binwidth = 0.5) +
  ggtitle("Replicated datasets of LYVE1 compared to the original data")
for (j in 1:4) {
  df <- data.frame(values=sample(draws$`ypred_1[2]`,(stan_data$N1+stan_data$N2+stan_data$N3)))
  pPostCheck_L <- pPostCheck_L + geom_histogram(data=df,aes(x=values),fill=j+1, alpha = 0.2, binwidth = 0.5)
}
pPostCheck_L

pPostCheck_R <- ggplot() + geom_histogram(data=data,aes(x=log(REG1B)),fill='white',color="black", binwidth = 0.5) +
  ggtitle("Replicated datasets of REG1B compared to the original data")
for (j in 1:4) {
  df <- data.frame(values=sample(draws$`ypred_1[3]`,(stan_data$N1+stan_data$N2+stan_data$N3)))
  pPostCheck_R <- pPostCheck_R + geom_histogram(data=df,aes(x=values),fill=j+1, alpha = 0.2, binwidth = 0.5)
}
pPostCheck_R

pPostCheck_T <- ggplot() + geom_histogram(data=data,aes(x=log(TFF1)),fill='white', color="black", binwidth = 0.5) +
  ggtitle("Replicated datasets of TFF1 compared to the original data")
for (j in 1:4) {
  df <- data.frame(values=sample(draws$`ypred_1[4]`,(stan_data$N1+stan_data$N2+stan_data$N3)))
  pPostCheck_T <- pPostCheck_T + geom_histogram(data=df,aes(x=values),fill=j+1, alpha = 0.2, binwidth = 0.5)
}
pPostCheck_T

