#Sensitivity analysis

library(rstan)
library(ggplot2)
data_path <- "Kurssit/Bayesian Data Analysis/Project/Debernardi et al 2020 data.csv" # Replace with working path
data <- read.csv(data_path)

# Split the data into test subject groups
group_1 <- data[data$diagnosis == 1,]
group_2 <- data[data$diagnosis == 2,]
group_3 <- data[data$diagnosis == 3,]

#Pooled

# Organize input data for the Stan model 
stan_data <- list(
  N1 = length(group_1[,1]),
  N2 = length(group_2[,1]),
  N3 = length(group_3[,1]),
  y1 = group_1[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y2 = group_2[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y3 = group_3[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  p = c(0,0,0,0,0,0,0,0)
)

#Parameters are given as p = [alpha, beta, alpha, beta, mu, sigma, alpha, beta]

#Scenarios
S <- matrix(c(1,1,0.5,1,0,20,1,1,
             2,2,1,2,0,1,2,2,
             0.01,0.01,0.01,0.01,10,1,0.01,0.01,
             100,100,100,100,1000,1000,100,100,
             100,0.01,100,0.01,-1000,1,100,0.01), ncol = 8, byrow = T)

p1 <- ggplot() + ggtitle('Creatinine')
p2 <- ggplot() + ggtitle('LYVE1')
p3 <- ggplot() + ggtitle('REG1B')
p4 <- ggplot() + ggtitle('TFF1')

for (i in 1:5) {
  # Fit the Stan model
  stan_data$p <- S[i,]
  pooled_fit <- stan(
    file = "Kurssit/Bayesian Data Analysis/Project/cancer-biomarker-analysis/pooled_sensitivity_analysis.stan",
    data = stan_data
  )
  monitor(pooled_fit)
  draws <- as.data.frame(pooled_fit)
  alpha <- mean(draws[,1])
  beta <- mean(draws[,2])
  values <- data.frame(x = seq(0,5,0.01), y=dgamma(seq(0,5,0.01), alpha, beta))
  p1 <- p1 + geom_line(aes(x, y), data = values, size = 1, col=i)
  mu <- mean(draws[,3])
  sigma <- mean(draws[,6])
  values <- data.frame(x = seq(-10,10,0.01), y=dnorm(seq(-10,10,0.01), mu, sigma))
  p2 <- p2 + geom_line(aes(x, y), data = values, size = 1, col=i)
  mu <- mean(draws[,4])
  sigma <- mean(draws[,7])
  values <- data.frame(x = seq(-5,15,0.01), y=dnorm(seq(-5,15,0.01), mu, sigma))
  p3 <- p3 + geom_line(aes(x, y), data = values, size = 1,col=i)
  mu <- mean(draws[,5])
  sigma <- mean(draws[,8])
  values <- data.frame(x = seq(-5,15,0.01), y=dnorm(seq(-5,15,0.01), mu, sigma))
  p4 <- p4 + geom_line(aes(x, y), data = values, size = 1, col=i)
}

p1
p2
p3
p4

# Hierarcical

# Organize input data for the Stan model 
stan_data <- list(
  N1 = length(group_1[,1]),
  N2 = length(group_2[,1]),
  N3 = length(group_3[,1]),
  y1 = group_1[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y2 = group_2[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  y3 = group_3[, c('creatinine', 'LYVE1', 'REG1B', 'TFF1')],
  p = c(0,0,0,0,0,0,0,0,0,0)
)

#Parameters are given as p = [alpha, beta, alpha, beta, alpha, beta, mu, sigma, alpha, beta]

#Scenarios
S <- matrix(c(1,1,1,1,1,1,0,20,1,1,
              2,2,2,2,2,2,0,1,2,2,
              0.01,0.01,0.01,0.01,0.01,0.01,10,1,0.01,0.01,
              100,100,100,100,100,100,1000,1000,100,100,
              100,0.01,100,0.01,100,0.01,-1000,1,100,0.01), ncol = 10, byrow = T)

p1 <- ggplot() + ggtitle('Creatinine')
p2 <- ggplot() + ggtitle('LYVE1')
p3 <- ggplot() + ggtitle('REG1B')
p4 <- ggplot() + ggtitle('TFF1')

for (i in 1:5) {
  # Fit the Stan model
  stan_data$p <- S[i,]
  hier_fit <- stan(
    file = "Kurssit/Bayesian Data Analysis/Project/cancer-biomarker-analysis/hierarchical_sensitivity_analysis.stan",
    data = stan_data
  )
  monitor(hier_fit)
  draws <- as.data.frame(hier_fit)
  
  for (j in seq(1,3)) {
    alpha <- mean(draws[,16+j])
    beta <- mean(draws[,28+j])
    values <- data.frame(x = seq(0,5,0.01), y=dgamma(seq(0,5,0.01), alpha, beta))
    p1 <- p1 + geom_line(aes(x, y), data = values, size = 1, col=j)
  }
  
  for (j in seq(1,3)) {
    mu <- mean(draws[,19+j])
    sigma <- mean(draws[,31+j])
    values <- data.frame(x = seq(-10,10,0.01), y=dnorm(seq(-10,10,0.01), mu, sigma))
    p2 <- p2 + geom_line(aes(x, y), data = values, size = 1,col=j)
  }
  
  for (j in seq(1,3)) {
    mu <- mean(draws[,22+j])
    sigma <- mean(draws[,34+j])
    values <- data.frame(x = seq(-10,10,0.01), y=dnorm(seq(-10,10,0.01), mu, sigma))
    p3 <- p3 + geom_line(aes(x, y), data = values, size = 1, col=j)
  }
  
  for (j in seq(1,3)) {
    mu <- mean(draws[,25+j])
    sigma <- mean(draws[,37+j])
    values <- data.frame(x = seq(-5,15,0.01), y=dnorm(seq(-5,15,0.01), mu, sigma))
    p4 <- p4 + geom_line(aes(x, y), data = values, size = 1, col=j)
  }

}

p1
p2
p3
p4

