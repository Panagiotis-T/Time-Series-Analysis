rm(list=ls())

## Predicting Lap Times in a Formula 1 Race
##Using the mathematical formulation for OLS and WLS. In the end, the relaxation method is used to find the optimal covariance matrix

##Import Plot Data
library(ggplot2)

##load data
load("D:/Documents/DTU/02417 Time Series Analysis/Assisgnment 1/DataA1.RData")
df <- AllDat
df = df[df$driverId == 9, ]



##OLS model


##Split the data into train [first 51 laps] and test [last 7 laps]

train = subset(AllDat, AllDat$lap_number<=51)
test = subset(AllDat, AllDat$lap_number>51)


coeff_mat <- data.frame(driverId=double(), theta0=double(), theta1=double())

##Build the matrices that comprise the design matrix

a_mat = matrix(0,dim(train)[1],max(train$driverId))

for(i in 1:dim(train)[1])  {
  idx = train[i,1]
  a_mat[i,idx] = 1
}


b_mat = matrix(train$lap_number)

c_mat = matrix(0,dim(train)[1])
for(i in 1:dim(train)[1])  {
  if (train[i,4]==TRUE) {
    c_mat[i] = 1
  }
}

d_mat = matrix(0,dim(train)[1])
for(i in 1:dim(train)[1])  {
  if (train[i,5]==1) {
    d_mat[i] = 1
  }
}

e_mat = matrix(0,dim(train)[1])

for(i in 1:dim(train)[1])  {
  if (train[i,2]==1) {
    e_mat[i] = 1
  }
}

##Combine the matrices / build the design matrix

X = cbind(a_mat,b_mat,c_mat,d_mat,e_mat)
Y = matrix(train$lap_time)

Y = train$lap_time

##Estimation of the parameters - theta 1 = 21 intercepts (for the 21 drivers), theta 22 = lap number, theta 23 = pit stop, theta 24 = pit stop lagged, theta 25 = first lap

##Calculate the parameters
theta <- solve(t(X)%*%X,t(X)%*%Y)

n <- length(Y)
p <- dim(X)[2]
err = Y-X%*%theta
sigma.hat <- sqrt(sum(err^2)/(n-p))

## Estimated parameters for driver 9

theta_9 <- theta[c(9,22,23,24,25)]
Y_9 <- subset(train, train$driverId==9)
Y9 <- t(t(Y_9$lap_time))
X9 <- X[,c(9,22,23,24,25)]
X9 = X[X[,9]==1,c(9,22,23,24,25)]

theta_9 = t(t(theta_9))
data9 <- subset(df,df$lap_number<=51)

##Variance of the parameters for driver 9
Var.theta_9 <- diag(sigma.hat^2*solve(t(X9)%*%X9))

##standard deviation of the parameters for driver 9
std.theta_9 <- sqrt(Var.theta_9)


## Build a new design matrix based on the test set

a_mat2 = matrix(0,dim(test)[1],max(train$driverId))

for(i in 1:dim(test)[1])  {
  idx = test[i,1]
  a_mat2[i,idx] = 1
}


b_mat2 = matrix(test$lap_number)

c_mat2 = matrix(0,dim(test)[1])
for(i in 1:dim(test)[1])  {
  if (test[i,4]==TRUE) {
    c_mat2[i] = 1
  }
}


d_mat2 = matrix(0,dim(test)[1])
for(i in 1:dim(test)[1])  {
  if (test[i,5]==1) {
    d_mat2[i] = 1
  }
}

e_mat2 = matrix(0,dim(test)[1])

for(i in 1:dim(test)[1])  {
  if (test[i,2]==1) {
    e_mat2[i] = 1
  }
}

## Predict lap time of the last 7 laps for driver 9

X_pre = cbind(a_mat2,b_mat2,c_mat2,d_mat2,e_mat2)
X9_pre = X_pre[X_pre[,9]==1,c(9,22,23,24,25)]
Y_pre = X9_pre%*%theta_9

## Measure of uncertainty /variance of the predictions

Var_pre9 <- diag(sigma.hat^2 * (matrix(1,7,7)+X9_pre%*%solve(t(X9)%*%X9)%*%t(X9_pre)))


## Calculate the 95% prediction interval

interval_pre9_test <- cbind(Y_pre,Y_pre) +cbind(matrix(-1,7,1) * t(t(qt(0.95,46)*sqrt(Var_pre9))),matrix(1,7,1)*t(t(qt(0.95,46)*sqrt(Var_pre9))))

Interval_pre9_train <- matrix(0,51,2)

Y9_train_pre <- X9%*%theta_9

for(i in 1:51)  {
  Un_root <- 1 + t(X9[i,])%*%solve(t(X9)%*%X9)%*%X9[i,]
  Interval_pre9_train[i,1] <- Y9_train_pre[i] - qt(0.95,46)*sigma.hat*sqrt(Un_root)
  Interval_pre9_train[i,2] <- Y9_train_pre[i] + qt(0.95,46)*sigma.hat*sqrt(Un_root)

}  

Interv_total <- as.data.frame(rbind(Interval_pre9_train,interval_pre9_test))


## Plot the fitted values, values for the observed lap times and 95% pred.interval

Y_fit <- as.data.frame(rbind(Y9_train_pre,Y_pre))
    
##ggplot

p1 <- ggplot()+
  geom_point(data=df, aes(x = lap_number, y= lap_time, color = "Observed \nvalues"), linetype="solid", size = 0.8)+
  geom_line(data=Y_fit, aes(x=as.integer(rownames(Y_fit)), y=V1, color="Fitted \nvalues"), linetype="solid", size = 0.7) +
  geom_line(data=Interv_total, aes(x=as.integer(rownames(Interv_total)), y=V1, color="Prediction \ninterval"), linetype="twodash", size = 1) +
  geom_line(data=Interv_total, aes(x=as.integer(rownames(Interv_total)), y=V2, color="Prediction \ninterval"), linetype="twodash", size = 1)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (12),hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica", size = (8)), 
                      axis.title = element_text(family = "Helvetica", size = (8), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "steelblue4", size = (8)),
                      legend.position = "bottom")


p1 + mynamestheme + labs( title= "Fitted and observed values with 95% prediction interval", y="Lap time (ms)", x = "Lap number", color = " ")
  


##WLS

## Build of the covariance matrix, assuming a correlation (rho) = 0.5

rho = rep(0.5, each=21)
LapsPerDriver <- c()
for(k in unique(AllDat$driverId)){
  LapsPerDriver <- c(LapsPerDriver,sum(AllDat$driverId==k))
}
LapsPerDriver_train <- pmin(LapsPerDriver,51)
library(magic)

##initialize Sigma matrix
Sigma <- matrix(0,ncol=0,nrow=0)
for(k in seq_along(LapsPerDriver_train)){
  N <- LapsPerDriver_train[k]
  O <- matrix(data = rep(0:(N-1),N),ncol=N)
  pow <- abs(t(O)-O)
  Block <- rho[k]^pow
  Sigma <- adiag(Sigma,Block)
}
##Covariance matrix
Sigma1 = Sigma

## Calculate the parameters

theta_w = solve(t(X)%*%solve(Sigma1)%*%X)%*%t(X)%*%solve(Sigma1)%*%Y

## Standard deviation of the parameters
sigma.hat_w = sqrt(1/(n-p)*t(Y-X%*%theta_w)%*%solve(Sigma1)%*%(Y-X%*%theta_w))


## Parameters for driver 9
theta_w_9 <- theta_w[c(9,22,23,24,25)]

## Sigma matrix for driver 9
Sigma9 <- matrix(0,ncol=0,nrow=0)
k=9
N <- LapsPerDriver_train[k]
O <- matrix(data = rep(0:(N-1),N),ncol=N)
pow <- abs(t(O)-O)
Sigma9 <- rho[k]^pow

## Variance of the parameters
Var.theta_9_w <- diag(as.numeric(sigma.hat_w)^2*solve(t(X9)%*%solve(Sigma9)%*%X9))

# Predict lap time based on new thetas for driver 9
Y_w_9 <- as.data.frame(X9%*%theta_w_9)

##Plot new lap time estimations

#ggplot
p2 <- ggplot()+
  geom_line(data=data9, aes(x = lap_number, y= lap_time, color = "Observed \nvalues"), linetype="solid", size = 0.8)+
  geom_line(data=Y_w_9, aes(x=as.integer(rownames(Y_w_9)), y=V1, color="Fitted \nvalues"), linetype="solid", size = 0.8)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (12),hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica", size = (8)), 
                      axis.title = element_text(family = "Helvetica", size = (8), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "steelblue4", size = (8)),
                      legend.position = "bottom")


p2 + mynamestheme + labs( title= "Fitted and observed values", y="Lap time (ms)", x = "Lap number", color = " ")


## Update covariance matrix by increasing the standard deviation 
## by a factor 1.5 for the first lap of each driver

gamma_pit<- 1.5
idx_train<- which(AllDat$lap_number<=51)
tmp1 <- AllDat$pitstop[idx_train]
tmp <- xor(AllDat$lap_number[idx_train]==1,tmp1)
Sigma2 <- Sigma1*(gamma_pit^outer(tmp,tmp,FUN="+"))


## Calculate parameters
theta_w_2 = solve(t(X)%*%solve(Sigma2)%*%X)%*%t(X)%*%solve(Sigma2)%*%Y
## standard deviation of the parameters
sigma.hat_w_2 = sqrt(1/(n-p)*t(Y-X%*%theta_w_2)%*%solve(Sigma2)%*%(Y-X%*%theta_w_2))


## WLS - Optimal covariance matrix

## Using the relaxation method to estimate the correlation and three factors 
##(gamma_first_lap, gamma_pit_stop, gamma_pit_lag)

rho_new = matrix(1,1,6)
gamma_first_lap <- matrix(1,1,6)
gamma_pit_stop <- matrix(1,1,6)
gamma_pit_lag <- matrix(1,1,6)

theta_w_4 = theta_w_2

##run the relaxation method for 5 iteration and see that it converges

for (i in 1:6) {
  ##compute residuals
  res = Y - X%*%theta_w_4
  ##rho is given by the correlation between the residuals for t and t-1
  rho_new[i]  = cor(res[-1],res[-n])
  ##set up Sigma
  Sigma_3 <- matrix(0,ncol=0,nrow=0)
  for(k in seq_along(LapsPerDriver_train)){
    N <- LapsPerDriver_train[k]
    O <- matrix(data = rep(0:(N-1),N),ncol=N)
    pow <- abs(t(O)-O)
    Block <- rho_new[i]^pow
    Sigma_3 <- adiag(Sigma_3,Block)
  }
  
  ##standard deviation of the residuals for laps with pit stop, first lap and pit stop lag
  sigma_pit <- sd(res[train$pitstop == TRUE])
  sigma_pit_lag <- sd(res[train$pitstop_lagged==1])
  sigma_first_lap <- sd(res[train$lap_number==1])
  ##standard deviation of the residuals for normal laps
  sigma_denom <- sd(res[train$lap_number !=1 & train$pitstop != TRUE & train$pitstop_lagged != 1])
  
  ## find gammas
  gamma_first_lap[i] <- sigma_first_lap/sigma_denom
  gamma_pit_stop[i] <- sigma_pit/sigma_denom
  gamma_pit_lag[i] <- sigma_pit_lag/sigma_denom
  
  
  idx_train<- which(AllDat$lap_number<=51)
  tmp1 <- AllDat$pitstop[idx_train]
  tmp2 <- AllDat$lap_number[idx_train]==1
  tmp3 <- AllDat$pitstop_lagged[idx_train]==1
  #covariance matrix
  Sigma_4 <- Sigma_3*(gamma_pit_stop[i]^outer(tmp1,tmp1,FUN="+"))*(gamma_first_lap[i]^outer(tmp2,tmp2,FUN="*"))*(gamma_pit_lag[i]^outer(tmp3,tmp3,FUN="+"))
  #estimate parameters
  theta_w_4 = solve(t(X)%*%solve(Sigma_4)%*%X)%*%t(X)%*%solve(Sigma_4)%*%Y
  
  #measure of uncertainty / standard deviation for the parameters
  sigma.hat_w__4 = as.numeric(sqrt(1/(n-p)*t(Y-X%*%theta_w_4)%*%solve(Sigma_4)%*%(Y-X%*%theta_w_4)))
  
}


##final parameters for driver 9
theta_9_final = theta_w_4[c(9,22,23,24,25)]

##last 7 laps prediction for driver 9
final_pred_9 = X9_pre%*%theta_9_final

##measure of uncertainty for final parameters

##variance
Var.theta_w_final <- diag(as.numeric(sigma.hat_w__4)^2*solve(t(X)%*%solve(Sigma_4)%*%X))

##standard deviation
std.theta_w_final <- sqrt(Var.theta_w_final)

##variance of the estimations
Var_final_pred9 <- diag((sigma.hat_w__4)^2 * (matrix(1,7,7)+X9_pre%*%solve(t(X9)%*%X9)%*%t(X9_pre)))


##find the 95% prediction interval

interval_pred9.final_test <- cbind(final_pred_9,final_pred_9) +cbind(matrix(-1,7,1) * t(t(qt(0.95,46)*sqrt(Var_final_pred9))),matrix(1,7,1)*t(t(qt(0.95,46)*sqrt(Var_final_pred9))))

Interval_pred9.final_train <- matrix(0,58,2)

##first 58 laps prediction for driver 9
Y9_train_pred_final <- X9%*%theta_9_final

##covariance matrix containing the values for driver 9
Sigma_final <- Sigma_4[431:481, 431:481]

##final design matrix for driver 9
X_9_tel <- rbind(X9, X9_pre)

Y9_train_test_pre <- rbind(Y9_train_pred_final, final_pred_9)



for(i in 1:58)  {
  un_root <- 1 + t(X_9_tel[i,])%*%solve(t(X9)%*%solve(Sigma_final)%*%X9)%*%X_9_tel[i,]
  Interval_pred9.final_train[i,1] <- Y9_train_test_pre[i] - qt(0.95,46)*sigma.hat_w__4*sqrt(un_root)
  Interval_pred9.final_train[i,2] <- Y9_train_test_pre[i] + qt(0.95,46)*sigma.hat_w__4*sqrt(un_root)
  
}

## fitted WLS
Y9_train_test_pre <- as.data.frame(Y9_train_test_pre)

Interv_total.final <- as.data.frame(Interval_pred9.final_train)

##fitted OLS
Y_ols <- as.data.frame(rbind(X9%*%theta_9,Y_pre))




#Plot fitted values from last model - 95% confidence interval, OLS fitted, observed values

#ggplot
p_f <- ggplot() + 
  geom_point(data=df, aes(x = lap_number, y= lap_time, color = "Observed \nlap times"), size = 0.8) + 
  geom_line(data=Interv_total.final, aes(x=as.integer(rownames(Interv_total.final)), y=V1, color="Prediction \ninterval"), linetype="twodash", size = 1) + 
  geom_line(data=Interv_total.final, aes(x=as.integer(rownames(Interv_total.final)), y=V2, color="Prediction \ninterval"), linetype="twodash", size = 1) +
  geom_line(data=Y9_train_test_pre, aes(x=as.integer(rownames(Y9_train_test_pre)), y=V1, color="WLS"), size = 0.7) +
  geom_line(data=Y_ols, aes(x=as.integer(rownames(Y_ols)), y=V1, color="OLS"), size = 0.6)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (9),hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica", size = (8)), 
                      axis.title = element_text(family = "Helvetica", size = (8), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "steelblue4", size = (8)),
                      legend.position = "bottom")

p_f + mynamestheme + labs( title= "Fitted values, OLS, observed lap times with 95% prediction interval", y="Lap time (ms)", x = "Lap number", color = " ")


