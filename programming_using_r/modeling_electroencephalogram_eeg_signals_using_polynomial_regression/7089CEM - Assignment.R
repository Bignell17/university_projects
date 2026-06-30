#Sets directory to where datasets are kept in my directory
setwd("C:\Users\Jordan\Desktop\Uni work\Intro to Stats for Data Science\7089 Assignment\7089CEM Assignment")

#Reading the EEG signals
x = read.table("x.csv", header = FALSE, sep = ",", dec = ".")
y = read.table("y.csv", header = FALSE, sep = ",", dec = ".")
time = read.table("time.csv", header = FALSE, sep = ",", dec = ".")

#renaming the column names
colnames(x) = paste0(rep("x",ncol(x)),1:ncol(x))
colnames(y) = "y"
colnames(time) = "time"

time = read.csv("time.csv", header = FALSE, skip = 1)
time = rbind(0, time)

#Defining EEG signals with and without time being used
df = cbind(x,y,time)
xy = cbind(x,y)

colnames(df)[which(names(df) == "V1")] <- "time"

#Libraries used
library(fitdistrplus)         #1.2
library(logspline)            #1.2
library(PerformanceAnalytics) #1.3
library(psych)                #1.3
library(ggExtra)              #3
library(ggplot2)              #3


##############################  Task 1  ##############################
############################## Task 1.1 ##############################
tsdf <- ts(df)    # Create a timeseries using all variables in df
tsdf              # Prints the time series
plot(tsdf)        # Plots the time series

############################## Task 1.2 ##############################
#defining each EEG signal
d1 <- density(df$x1)
d2 <- density(df$x2)
d3 <- density(df$x3)
d4 <- density(df$x4)
d5 <- density(df$y)

#Using Cullen & Frey graph to identify distribution type and other measurements
descdist(df$y, discrete = FALSE)
descdist(df$x1, discrete = FALSE)
descdist(df$x2, discrete = FALSE)
descdist(df$x3, discrete = FALSE)
descdist(df$x4, discrete = FALSE)

#Plotting and showing distribution for each EEG signal
windows(8,8)
ylim1=c(0,0.5)
plot(d1, col="blue", ylim=ylim1, xlab="EEG Signals", lwd=2, main = "Distribution for each EEG Signal")
lines(d2, col="red", lwd=2)
lines(d3, col="green",lwd=2)
lines(d4, col="yellow",lwd=2)
lines(d5, col="pink", lwd=2)
legend("topright",
       legend=c("x1", "x2", "x3", "x4", "y"),
       col=c("blue", "red", "green", "yellow", "pink"),
       text.col= "Black",
       lwd=2)

############################## Task 1.3 ##############################

#Plotting using correlation and scatterplot to show relations and dependencies
chart.Correlation(xy,histogram = TRUE)

library(psych)

describe(df)
summary(df)


##############################  Task 2  ##############################

#Creating and following the correct structure for each model
ones <- matrix(1,201,1)
model1 <- cbind(ones, df$x4, df$x1^2, df$x1^3, df$x3^4)
model2 <- cbind(ones, df$x3^3, df$x3^4)
model3 <- cbind(ones, df$x2, df$x1^3, df$x3^4)
model4 <- cbind(ones, df$x4, (df$x1^3), (df$x3^4))
model5 <- cbind(ones, df$x4, df$x1^2, (df$x1^3), (df$x3^4), (df$x1^4))
  
############################## Task 2.1 ##############################

#Estimating the theta hat for each EEG signal
theta_hat1 <- solve(t(model1) %*% model1) %*% t(model1) %*% df$y

#Printing the result of each theta hat
{
  print(paste("Theta Bias: ", theta_hat1[1]))
  print(paste("Theta 1: ", theta_hat1[2]))
  print(paste("Theta 2: ", theta_hat1[3]))
  print(paste("Theta 3: ", theta_hat1[4]))
  print(paste("Theta 4: ", theta_hat1[5]))
}

theta_hat2 <- solve(t(model2) %*% model2) %*% t(model2) %*% df$y
{
  print(paste("Theta Bias: ",theta_hat2[1]))
  print(paste("Theta 1: ",theta_hat2[2]))
  print(paste("Theta 2: ",theta_hat2[3]))
}

theta_hat3 <- solve(t(model3) %*% model3) %*% t(model3) %*% df$y
{
  print(paste("Theta Bias: ",theta_hat3[1]))
  print(paste("Theta 1: ",theta_hat3[2]))
  print(paste("Theta 2: ",theta_hat3[3]))
  print(paste("Theta 3: ",theta_hat3[4]))
}

theta_hat4 <- solve(t(model4) %*% model4) %*% t(model4) %*% df$y
{
  print(paste("Theta Bias: ",theta_hat4[1]))
  print(paste("Theta 1: ",theta_hat4[2]))
  print(paste("Theta 2: ",theta_hat4[3]))
  print(paste("Theta 3: ",theta_hat4[4]))
}

theta_hat5 <- solve(t(model5) %*% model5) %*% t(model5) %*% df$y
{
  print(paste("Theta Bias: ",theta_hat5[1]))
  print(paste("Theta 1: ",theta_hat5[2]))
  print(paste("Theta 2: ",theta_hat5[3]))
  print(paste("Theta 3: ",theta_hat5[4]))
  print(paste("Theta 4: ",theta_hat5[5]))
  print(paste("Theta 5: ",theta_hat5[6]))
}


############################## Task 2.2 ##############################

#Defining the Y hat before calculating the RSS models
yhat1 <- model1 %*% theta_hat1
yhat2 <- model2 %*% theta_hat2
yhat3 <- model3 %*% theta_hat3
yhat4 <- model4 %*% theta_hat4
yhat5 <- model5 %*% theta_hat5

#Calculating the RSS model
RSS_model1 <- sum((df$y - yhat1)^2)
RSS_model2 <- sum((df$y - yhat2)^2)
RSS_model3 <- sum((df$y - yhat3)^2)
RSS_model4 <- sum((df$y - yhat4)^2)
RSS_model5 <- sum((df$y - yhat5)^2)

#printing the RSS models
{
  print(paste("RSS Model 1:", RSS_model1))
  print(paste("RSS Model 2:", RSS_model2))
  print(paste("RSS Model 3:", RSS_model3))
  print(paste("RSS Model 4:", RSS_model4))
  print(paste("RSS Model 5:", RSS_model5))
}


############################## Task 2.3 ##############################

#Defining n as the number of data samples
n <- 201

#Computing sigma
sigma1 <- (RSS_model1/(n-1))
sigma2 <- (RSS_model2/(n-1))
sigma3 <- (RSS_model3/(n-1))
sigma4 <- (RSS_model4/(n-1))
sigma5 <- (RSS_model5/(n-1))

#Computing the log-likelihood for each model
logp_model1 <- -(201/2) * log(2*pi) - (201/2) * log(sigma1) - (1/(2*sigma1)) * RSS_model1
logp_model2 <- -(201/2) * log(2*pi) - (201/2) * log(sigma2) - (1/(2*sigma2)) * RSS_model2
logp_model3 <- -(201/2) * log(2*pi) - (201/2) * log(sigma3) - (1/(2*sigma3)) * RSS_model3
logp_model4 <- -(201/2) * log(2*pi) - (201/2) * log(sigma4) - (1/(2*sigma4)) * RSS_model4
logp_model5 <- -(201/2) * log(2*pi) - (201/2) * log(sigma5) - (1/(2*sigma5)) * RSS_model5

#Printing each log-likelihood model
{
  print(paste("Log-likelihood Model 1:", logp_model1))
  print(paste("Log-likelihood Model 2:", logp_model2))
  print(paste("Log-likelihood Model 3:", logp_model3))
  print(paste("Log-likelihood Model 4:", logp_model4))
  print(paste("Log-likelihood Model 5:", logp_model5))
}


############################## Task 2.4 ##############################

#defining the estimated parameters for each model
{
  k1<-5 #Length of model 1
  k2<-3 #Length of model 2
  k3<-4 #Length of model 3
  k4<-4 #Length of model 4
  k5<-6 #Length of model 5
}
#Calculating the AIC/BIC models
{
  aicmodel1 <- (2*k1-2*logp_model1)
  aicmodel2 <- (2*k2-2*logp_model2)
  aicmodel3 <- (2*k3-2*logp_model3)
  aicmodel4 <- (2*k4-2*logp_model4)
  aicmodel5 <- (2*k5-2*logp_model5)

  bicmodel1 <- (k1*log(n) - 2*logp_model1)
  bicmodel2 <- (k2*log(n) - 2*logp_model2)
  bicmodel3 <- (k3*log(n) - 2*logp_model3)
  bicmodel4 <- (k4*log(n) - 2*logp_model4)
  bicmodel5 <- (k5*log(n) - 2*logp_model5)
}
#Printing the AIC/BIC models in a table format
{
  AICBIC = matrix(cbind(aicmodel1,aicmodel2,aicmodel3,aicmodel4,aicmodel5,
                              bicmodel1,bicmodel2,bicmodel3,bicmodel4,bicmodel5),
                        nrow=5,
                        ncol=2)
  colnames(AICBIC) <- c("AIC", "BIC") 
  rownames(AICBIC) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5" )
  AICBIC
}

############################## Task 2.5 ##############################
#Computing the residuals for each model
{
  resid1 <- (df$y - yhat1)
  resid2 <- (df$y - yhat2)
  resid3 <- (df$y - yhat3)
  resid4 <- (df$y - yhat4)
  resid5 <- (df$y - yhat5)
}

#Plotting each QQ-Plot model

windows(10,10)
{
  qqnorm(resid1, pch=1, frame =FALSE)
  qqline(resid1, col="blue", lwd=2)
}
{
  qqnorm(resid2, pch=1, frame =FALSE)
  qqline(resid2, col="red", lwd=2)
}
{
  qqnorm(resid3, pch=1, frame =FALSE)
  qqline(resid3, col="green", lwd=2)
}
{
  qqnorm(resid4, pch=1, frame =FALSE)
  qqline(resid4, col="yellow", lwd=2)
}
{
  qqnorm(resid5, pch=1, frame =FALSE)
  qqline(resid5, col="purple", lwd=2)
}


############################## Task 2.6 ##############################
#Plotting the AIC/BIC models from task 2.4 and printing the best model for each AIC/BIC
{
  AICBIC = matrix(cbind(aicmodel1,aicmodel2,aicmodel3,aicmodel4,aicmodel5,
                     bicmodel1,bicmodel2,bicmodel3,bicmodel4,bicmodel5),
                     nrow=5,
                     ncol=2)
  colnames(AICBIC) <- c("AIC", "BIC") 
  rownames(AICBIC) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5" )
  
  AIC = cbind(aicmodel1,aicmodel2,aicmodel3,aicmodel4,aicmodel5)
  BIC = cbind(bicmodel1,bicmodel2,bicmodel3,bicmodel4,bicmodel5)
  bestAICmod <- apply(AIC,1,min)
  bestBICmod <- apply(BIC,1,min)
  print(AICBIC)
  print(paste("Best AIC model is model 3:", bestAICmod))
  print(paste("Best BIC model is model 3:", bestBICmod))
}

############################## Task 2.7 ##############################

#Splitting the EEG signals X and Y into a test and train with 70/30 split
rows=nrow(xy)
index=0.70*rows
train = df[1:index, ]
test = df[-(1:index),]

train
test


##############################  2.7.1  ###############################

#Creating the model and using train data
onesTrain <- matrix(1,140,1)
xTrain <- cbind(onesTrain, train$x2, train$x1^3, train$x3^4)

#Finding the theta hat from train data
ThetaHatTrain <- solve(t(xTrain) %*% xTrain) %*% t(xTrain) %*% train$y
{
  print(paste("Theta Bias: ",ThetaHatTrain[1]))
  print(paste("Theta 1: ",ThetaHatTrain[2]))
  print(paste("Theta 2: ",ThetaHatTrain[3]))
  print(paste("Theta 3: ",ThetaHatTrain[4]))
}
##############################  2.7.2  ###############################

#Creating the model and using test data
onesTest = matrix(1,61,1)
xTest <- cbind(onesTest, test$x2, test$x1^3, test$x3^4)

#Finding the y hat from test data
yHatTest <- xTest %*% ThetaHatTrain

#printing the y hat test
yHatTest

##############################  2.7.3  ###############################

#Defining RSS for test data
RSSTest = sum((test$y - yHatTest)^2)

#Defining Sigma
sigmaTest = RSSTest/( 61 - 1 )

#Defining Covariance 
cov_thetaHat = sigmaTest * (solve(t(xTest) %*% xTest))
var_y_hat = matrix(0 , 61 , 1)

#Computing the Confidence interval
for( i in 1:61)
{
  X_i = matrix(xTest[i,] , 1 , 4 ) # X[i,] creates a vector. Convert it to matrix
  var_y_hat[i,1] = X_i %*% cov_thetaHat %*% t(X_i) # same as sigma_2 * ( X_i %*% ( solve(t(X) %*% X)  ) %*% t(X_i) )
}

CI = 1.96 * sqrt(var_y_hat) # We use 1.96 because we compute the model prediction at 95% CI

#Defining the upper and lower confidence intervals
upperCI <- yHatTest+CI
lowerCI <- yHatTest-CI

#plotting the error bars
windows(20,10)
plot(test$time, yHatTest, ,ylim = c(-2,3))
segments(test$time, lowerCI, test$time, upperCI)
arrows(test$time, lowerCI, test$time, upperCI, length=0.05, angle=90, code=3)
lines(test$time, yHatTest,type = 'o',col = 'purple')



##############################  Task 3  ##############################
#Calculating the Largest absoulte value from theta hat 3
#Theta hat 3 is regarded 'best' model  
largeabsval <- matrix(cbind(theta_hat3))
lav1 <- largeabsval[1,]
lav2 <- largeabsval[2,]
lav3 <- cbind(lav1, lav2)
lav3

#Creating sample number 
sample_number = 5000

#Creating prior uniform distribution of 2 largest absolute values (theta Bias and theta 1)
priorunibias <- runif(sample_number, lav1 * (0.5), lav1 * (1.5))
prioruni1 <- runif(sample_number, lav2 * (0.5), lav2 * (1.5)) #min/max of numbers being included - range

#Creating new ABC theta hat
ABC_theta_hat = cbind(priorunibias, prioruni1, theta_hat3[3], theta_hat3[4])

#Creating variables for both accepted/rejected thetas in following computation
accepted_thetas = matrix(0, 0, 2)
rejected_thetas = matrix(0, 0, 2)

#Computing the distance against the epsilon
#Puts data points into categories of accepted/rejected thetas dependable from calculation
for (i in 1:sample_number){
  if(sum(df$y - model3 %*% ABC_theta_hat[i,])^2 < (1.5 * RSS_model3)){ #Calculating the distance against the epsilon
    accepted_thetas = rbind(accepted_thetas,c(ABC_theta_hat[i,1], ABC_theta_hat[i,2])) #Accepted/checked thetas
  }
  else{
    rejected_thetas = rbind(rejected_thetas,c(ABC_theta_hat[i,1], ABC_theta_hat[i,2])) #Rejected thetas
  }
}

############## Plotting joint and marginal posterior distribution ##############

ABCgraph = ggExtra::ggMarginal(ggplot2::ggplot(data.frame(accepted_thetas), ggplot2::aes(x=priorunibias, y=prioruni1))+
                           ggplot2::geom_point(colour = "green")+
                           ggplot2::geom_point(data = data.frame(rejected_thetas), colour = "red", alpha = 0.2),
                         type = "densigram",
                         fill = "green")
ABCgraph


##########################  Contour plot  ############################

ABCcontour = ggplot(data.frame(accepted_thetas), aes(x=priorunibias, y=prioruni1))+
   stat_density_2d_filled(contour = TRUE, geom = "polygon", colour="white")

ABCcontour


