#Plot histograms by groups

data_path <- "Kurssit/Bayesian Data Analysis/Project/Debernardi et al 2020 data.csv" #'~/your/path/to/dataset' # Replace with working path
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
