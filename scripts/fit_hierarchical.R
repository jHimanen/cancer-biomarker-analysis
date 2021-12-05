##################################################
# This is an R script used to fit the
# hierarchical Stan model and compute convergence, 
# model evaluation, and posterior predictive
# diagnostics for the "Urinary
# biomarkers for pancreatic cancer" dataset.
##################################################

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
hier_fit <- stan(
  file = '.../cancer-biomarker-analysis/models/hierarchical.stan', # REPLACE WITH WORKING PATH
  data = stan_data,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

# Store convergence diagnostics
results <- monitor(hier_fit)
selected_res <- results[c(1:40), c('mean', 'sd', 'n_eff', 'Rhat', 'Q5', 'Q50', 'Q95')]
res_df <- as.data.frame(selected_res)
write.csv(res_df,'.../cancer-biomarker-analysis/diagnostic_data/hier_res.csv') # REPLACE WITH WORKING PATH

# Extract log-likelihoods for LOO evaluation
log_liks <- list(
  extract_log_lik(hier_fit, parameter_name = 'log_lik_1', merge_chains = FALSE),
  extract_log_lik(hier_fit, parameter_name = 'log_lik_2', merge_chains = FALSE),
  extract_log_lik(hier_fit, parameter_name = 'log_lik_3', merge_chains = FALSE)
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

hier_eval <- as.data.frame(diagnostics)
write.csv(hier_eval, '.../cancer-biomarker-analysis/diagnostic_data/hier_eval.csv') # REPLACE WITH WORKING PATH

# Visualize posterior distributions

draws <- as.data.frame(hier_fit)

pCreatinine <- ggplot() + ggtitle('Posterior distribution of creatinine')
x <- seq(0,5,0.05)
for (i in seq(1,8000,100)) {
  df1 <- data.frame(x=x,y=dgamma(x,draws$`alpha[1]`[i],draws$`beta[1]`[i]))
  df2 <- data.frame(x=x,y=dgamma(x,draws$`alpha[2]`[i],draws$`beta[2]`[i]))
  df3 <- data.frame(x=x,y=dgamma(x,draws$`alpha[3]`[i],draws$`beta[3]`[i]))
  pCreatinine <- pCreatinine + 
    geom_line(data=df1,aes(x,y), color="lightgreen",alpha=0.4) +    
    geom_line(data=df2,aes(x,y), color="skyblue",alpha=0.4) +
    geom_line(data=df3,aes(x,y), color="tomato",alpha=0.4)
}
pCreatinine

pLYVE1 <- ggplot() + ggtitle('Posterior distribution of LYVE1')
x <- seq(-15,15,0.1)
for (i in seq(1,8000,100)) {
  df1 <- data.frame(x=x,y=dnorm(x,draws$`mu[1,1]`[i],draws$`sigma[1,1]`[i]))
  df2 <- data.frame(x=x,y=dnorm(x,draws$`mu[1,2]`[i],draws$`sigma[1,2]`[i]))
  df3 <- data.frame(x=x,y=dnorm(x,draws$`mu[1,3]`[i],draws$`sigma[1,3]`[i]))
  pLYVE1 <- pLYVE1 + 
    geom_line(data=df1,aes(x,y),col='lightgreen',alpha=0.4) +    
    geom_line(data=df2,aes(x,y),col='skyblue',alpha=0.4) +
    geom_line(data=df3,aes(x,y),col='tomato',alpha=0.4)
}
pLYVE1

pREG1B <- ggplot() + ggtitle('Posterior distribution of REG1B')
x <- seq(-10,15,0.1)
for (i in seq(1,8000,100)) {
  df1 <- data.frame(x=x,y=dnorm(x,draws$`mu[2,1]`[i],draws$`sigma[2,1]`[i]))
  df2 <- data.frame(x=x,y=dnorm(x,draws$`mu[2,2]`[i],draws$`sigma[2,2]`[i]))
  df3 <- data.frame(x=x,y=dnorm(x,draws$`mu[2,3]`[i],draws$`sigma[2,3]`[i]))
  pREG1B <- pREG1B + 
    geom_line(data=df1,aes(x,y),col='lightgreen',alpha=0.4) +    
    geom_line(data=df2,aes(x,y),col='skyblue',alpha=0.4) +
    geom_line(data=df3,aes(x,y),col='tomato',alpha=0.4)
}
pREG1B

pTFF1 <- ggplot() + ggtitle('Posterior distribution of TFF1')
x <- seq(-3,15,0.1)
for (i in seq(1,8000,100)) {
  df1 <- data.frame(x=x,y=dnorm(x,draws$`mu[3,1]`[i],draws$`sigma[3,1]`[i]))
  df2 <- data.frame(x=x,y=dnorm(x,draws$`mu[3,2]`[i],draws$`sigma[3,2]`[i]))
  df3 <- data.frame(x=x,y=dnorm(x,draws$`mu[3,3]`[i],draws$`sigma[3,3]`[i]))
  pTFF1 <- pTFF1 + 
    geom_line(data=df1,aes(x,y),col='lightgreen',alpha=0.4) +    
    geom_line(data=df2,aes(x,y),col='skyblue',alpha=0.4) +
    geom_line(data=df3,aes(x,y),col='tomato',alpha=0.4)
}
pTFF1

# Posterior predictive check

# Visualizing posterior predictive distributions
stan_hist(hier_fit, pars = c('ypred_1[1]','ypred_1[2]','ypred_1[3]','ypred_1[4]'))
stan_hist(hier_fit, pars = c('ypred_2[1]','ypred_2[2]','ypred_2[3]','ypred_2[4]'))
stan_hist(hier_fit, pars = c('ypred_3[1]','ypred_3[2]','ypred_3[3]','ypred_3[4]'))

plotPostCheck <- function(dataCol, ypred, N, bin, title, xlabel){
  df<- data.frame(values=dataCol)
  p <- ggplot() + 
    geom_histogram(data=df,aes(x=values),fill='white',color="black", binwidth = bin) +
    ggtitle(title) + xlab(xlabel)
  for (j in 1:4) {
    df <- data.frame(values=sample(ypred,N))
    p <- p + geom_histogram(data=df,aes(x=values),fill=j+1, alpha = 0.2, binwidth = bin)
  }
  p
}

#Creatinine
pPostCheck_C1 <- plotPostCheck(data[data$diagnosis==1,]$creatinine,draws$`ypred_1[1]`,stan_data$N1,0.1,
                              "Replicated datasets of creatinine compared to the original data for diagnosis 1",
                              "creatinine")
pPostCheck_C1
pPostCheck_C2 <- plotPostCheck(data[data$diagnosis==2,]$creatinine,draws$`ypred_2[1]`,stan_data$N2,0.1,
                               "Replicated datasets of creatinine compared to the original data for diagnosis 2 ",
                               "creatinine")
pPostCheck_C2
pPostCheck_C3 <- plotPostCheck(data[data$diagnosis==3,]$creatinine,draws$`ypred_3[1]`,stan_data$N3,0.1,
                               "Replicated datasets of creatinine compared to the original data for diagnosis 3",
                               "creatinine")
pPostCheck_C3

#LYVE1
pPostCheck_L1 <- plotPostCheck(log(data[data$diagnosis==1,]$LYVE1),draws$`ypred_1[2]`,stan_data$N1,0.5,
                               "Replicated datasets of LYVE1 compared to the original data for diagnosis 1",
                               "log(LYVE1)")
pPostCheck_L1
pPostCheck_L2 <- plotPostCheck(log(data[data$diagnosis==2,]$LYVE1),draws$`ypred_2[2]`,stan_data$N2,0.5,
                               "Replicated datasets of LYVE1 compared to the original data for diagnosis 2",
                               "log(LYVE1)")
pPostCheck_L2
pPostCheck_L3 <- plotPostCheck(log(data[data$diagnosis==3,]$LYVE1),draws$`ypred_3[2]`,stan_data$N3,0.5,
                               "Replicated datasets of LYVE1 compared to the original data for diagnosis 3",
                               "log(LYVE1)")
pPostCheck_L3

#REG1B
pPostCheck_R1 <- plotPostCheck(log(data[data$diagnosis==1,]$REG1B),draws$`ypred_1[3]`,stan_data$N1,0.5,
                               "Replicated datasets of REG1B compared to the original data for diagnosis 1",
                               "log(REG1B)")
pPostCheck_R1
pPostCheck_R2 <- plotPostCheck(log(data[data$diagnosis==2,]$REG1B),draws$`ypred_2[3]`,stan_data$N2,0.5,
                               "Replicated datasets of REG1B compared to the original data for diagnosis 2",
                               "log(REG1B)")
pPostCheck_R2
pPostCheck_R3 <- plotPostCheck(log(data[data$diagnosis==3,]$REG1B),draws$`ypred_3[3]`,stan_data$N2,0.5,
                               "Replicated datasets of REG1B compared to the original data for diagnosis 3",
                               "log(REG1B)")
pPostCheck_R3

#TFF1
pPostCheck_T1 <- plotPostCheck(log(data[data$diagnosis==1,]$TFF1),draws$`ypred_1[4]`,stan_data$N1,1,
                               "Replicated datasets of TFF1 compared to the original data for diagnosis 1",
                               "log(TFF1)")
pPostCheck_T1
pPostCheck_T2 <- plotPostCheck(log(data[data$diagnosis==2,]$TFF1),draws$`ypred_2[4]`,stan_data$N2,1,
                               "Replicated datasets of TFF1 compared to the original data for diagnosis 2",
                               "log(TFF1)")
pPostCheck_T2
pPostCheck_T3 <- plotPostCheck(log(data[data$diagnosis==3,]$TFF1),draws$`ypred_3[4]`,stan_data$N3,1,
                               "Replicated datasets of TFF1 compared to the original data for diagnosis 3",
                               "log(TFF1)")
pPostCheck_T3

