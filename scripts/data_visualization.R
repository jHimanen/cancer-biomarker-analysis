# Plot histograms by groups

library(ggplot2)

data_path <- '.../Debernardi et al 2020 data.csv' # REPLACE WITH WORKING PATH
data <- read.csv(data_path)

#Creatinine
ggplot(data,aes(x=creatinine)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 0.2) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 0.2) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 0.2)

#LYVE1
ggplot(data,aes(x=LYVE1)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 1)

data$LYVE1 <- log(data$LYVE1)
ggplot(data,aes(x=LYVE1)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 0.5) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 0.5) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 0.5)

#REG1B
ggplot(data,aes(x=REG1B)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 50) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 50) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 50)  

data$REG1B <- log(data$REG1B)
ggplot(data,aes(x=REG1B)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 1) 

#TFF1
ggplot(data,aes(x=TFF1)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 500) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 500) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 500)

data$TFF1 <- log(data$TFF1)
ggplot(data,aes(x=TFF1)) + 
  geom_histogram(data=subset(data,diagnosis==1),fill = "green", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==2),fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data=subset(data,diagnosis==3),fill = "red", alpha = 0.2, binwidth = 1)

#Plot prior distrbutions

x <- seq(0,5,0.01)
y1 <- dgamma(x,1,1)
y2 <- dgamma(x,0.5,1)
df <- data.frame(x,y1,y2)
ggplot(data=df) + geom_line(aes(x,y1, colour ="Gamma(1,1)")) + geom_line(aes(x,y2,colour="Gamma(0.5,1)")) +
  ylab('y') +
  scale_color_manual(name = "Priors", values = c("Gamma(1,1)" = "darkblue", "Gamma(0.5,1)" = "red"))+
  theme(legend.position=c(0.75,0.85))


