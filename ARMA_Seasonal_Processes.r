
#Panagiotis Tsagkaroulis

rm(list=ls())


library(tseries)

##Check stationarity of the process
zroots <- polyroot(c(1, 1, 1, 1))


#Simulate 10 realisations with 200 observations from the process and plot them
#install.packages('tidyr', repos='http://cran.rstudio.com/') 
sim1 <- cbind("Time"=1:200, replicate(10,arima.sim(model = list(ma=c(1,1,1), order=c(0,0,3)), sd = 0.1, n = 200)))
colnames(sim1) <- c('Time','1','2','3','4','5','6','7','8','9','10')
dfsim1 = tidyr::pivot_longer(as.data.frame(sim1), cols=2:11)
library(ggplot2)


p <- ggplot(dfsim1, aes(y=value, color=name, x=Time)) + geom_line()
p + theme_light()+theme(axis.title=element_text(size=10), legend.position = "none",
                        axis.title.x = element_blank())


#Estimate the ACF for each realisation and plot those
acfsim1 <- cbind("Lag"=0:24, apply(sim1[1:25,2:11], 2, function(x){acf(x, 25, plot=FALSE)$acf}))
dfacfsim1 = tidyr::pivot_longer(as.data.frame(acfsim1), cols=2:11)
#the following plots in classical ACF form with spikes
ggplot(dfsim1, aes(y=value, color=name, x=Lag)) + geom_segment(data=dfacfsim1, aes(xend=Lag, yend=0))
#the following plots by joining the points 
#ggplot(dfsim1, aes(y=value, color=name, x=Lag)) + geom_segment(data=dfacfsim1)
p1 <- ggplot(dfacfsim1, aes(y=value, color=name, x=Lag)) + geom_line()
p1 + theme_light()+theme(axis.title=element_text(size=10), legend.position = "none",
                        axis.title.x = element_blank())+ labs(y="ACF")+
  geom_text(aes(label=name))



#Repeat for the PACF of the same realisations
pacfsim1 <- cbind("Lag"=1:24, apply(sim1[1:25,2:11], 2, function(x){acf(x, 25, plot=FALSE, type=c("partial"))$acf}))
dfpacfsim1 = tidyr::pivot_longer(as.data.frame(pacfsim1), cols=2:11)

#the following plots in classical ACF form with spikes
ggplot(dfsim1, aes(y=value, color=name, x=Lag)) + geom_segment(data=dfpacfsim1, aes(xend=Lag, yend=0))
#the following plots by joining the points 
#ggplot(dfsim1, aes(y=value, color=name, x=Lag)) + geom_segment(data=dfacfsim1)
p2 <- ggplot(dfpacfsim1, aes(y=value, color=name, x=Lag)) + geom_line()
p2 + theme_light()+theme(axis.title=element_text(size=10), legend.position = "none",
                         axis.title.x = element_blank())+ labs(y="PACF")+
  geom_text(aes(label=name))



#Calculate the variance of each of the realisations
Variance <-apply(sim1[,2:11],2,var)
print(Variance)
print(mean(Variance))

##Predicting temperatures in a district heating network
#creating the data frame
year <- c(rep(2016,4),rep(2017,11))
month <- c(9:12,1:11)
temp <- c( 46.6, 49.5, 60.3, 59.2, 59.5, 61.9, 59.7, 60.1, 57.8, 49.7, 49.7, 50.1, 48.6, 54.5, 62.3)
df <- data.frame(year,month,temp)

#mu = 55
#transformation
X.t <- df$temp-55
X.len <- length(X.t)

#prediction 
X.hat2017M12 <- 0.5*X.t[X.len] - 0.3*X.t[X.len-1] + 0.9*X.t[X.len-11] - 0.45*X.t[X.len-12] +0.27*X.t[X.len-13]
X.hat2018M1 <- 0.5*X.hat2017M12 - 0.3*X.t[X.len] + 0.9*X.t[X.len-10] - 0.45*X.t[X.len-11] +0.27*X.t[X.len-12]

#backtransformation
Y.hat2017M12 <- X.hat2017M12 + 55
Y.hat2018M1 <- X.hat2018M1 + 55

#variance of prediction error (et+1 and et+2)
sigma.epsilon<-0.5
V.epsilon2017M12<-sigma.epsilon^2
V.epsilon2018M1<-(1+(0.5)^2)*(sigma.epsilon^2)

#95% confidence interval
inter<-qnorm(0.975)


Y2017M12.low<-Y.hat2017M12-inter*sqrt(V.epsilon2017M12)
Y2017M12.high<-Y.hat2017M12+inter*sqrt(V.epsilon2017M12)
Y2018M1.low<-Y.hat2018M1-inter*sqrt(V.epsilon2018M1)
Y2018M1.high<-Y.hat2018M1+inter*sqrt(V.epsilon2018M1)

#Plot


prediction.interval.low<-c(Y2017M12.low,Y2018M1.low)
prediction.interval.high<-c(Y2017M12.high,Y2018M1.high)


##save figure
#jpeg("pred_ci.jpeg", units="in", width=5, height=4, res=600)

time=c(1:15)
plot(time,df$temp, xlim=c(0,18),ylim=c(min(df$temp)-5,max(df$temp)+5),main=expression('95% Confidence interval',hjust = 0.5),
     xlab='Time', ylab='Temperature')
points(16,Y.hat2017M12,col = "red") 
points(17,Y.hat2018M1,col = "red")
points(16:17,prediction.interval.low,  col=2,pch=95, lty=c(1,2,2))
points(16:17,prediction.interval.high, col=2,pch=95, lty=c(1,2,2))
lines(time,df$temp, pch=16)
par(xpd=TRUE)
legend(x = "bottomright", inset = 0.01, legend = c("prediction", "conf.int"), pch = c(1, 95), col = c(2,2), cex=0.8)                

#dev.off()

##Simulating seasonal processes

#arima.sim function
arima.sim2 <- function (model, n, rand.gen = rnorm, innov = rand.gen(n, ...), 
          n.start = NA, start.innov = rand.gen(n.start, ...), ...) 
{
  if (!is.list(model)) 
    stop("'model' must be list")
  if (n <= 0L) 
    stop("'n' must be strictly positive")
  p <- length(model$ar)
  if (p) {
    minroots <- min(Mod(polyroot(c(1, -model$ar))))
    if (minroots <= 1) 
      warning ("'ar' part of model is not stationary")
  }
  q <- length(model$ma)
  if (is.na(n.start)) 
    n.start <- p + q + ifelse(p > 0, 
                              0)
  if (n.start < p + q) 
    stop("burn-in 'n.start' must be as long as 'ar + ma'")
  d <- 0
  if (!is.null(ord <- model$order)) {
    if (length(ord) != 3L) 
      stop("'model$order' must be of length 3")
    if (p != ord[1L]) 
      stop("inconsistent specification of 'ar' order")
    if (q != ord[3L]) 
      stop("inconsistent specification of 'ma' order")
    d <- ord[2L]
    if (d != round(d) || d < 0) 
      stop("number of differences must be a positive integer")
  }
  if (!missing(start.innov) && length(start.innov) < n.start) 
    stop(sprintf(ngettext(n.start, "'start.innov' is too short: need %d point", 
                          "'start.innov' is too short: need %d points"), 
                 n.start), domain = NA)
  x <- ts(c(start.innov[seq_len(n.start)], innov[1L:n]), start = 1 - 
            n.start)
  if (length(model$ma)) {
    x <- filter(x, c(1, model$ma), sides = 1L)
    x[seq_along(model$ma)] <- 0
  }
  if (length(model$ar)) 
    x <- filter(x, model$ar, method = "recursive")
  if (n.start > 0) 
    x <- x[-(seq_len(n.start))]
  if (d > 0) 
    x <- diffinv(x, differences = d)
  as.ts(x)
}




#ARIMA(1,0,0)x(0,0,0)12 phi1=-0.85 in reality it is an AR(1) model


##windows(5, 5)
sim3.1 <- arima.sim(model = list(ar=0.85, order=c(1,0,0)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.1,main=expression(''),ylab=' ',xlab=' ')
acf1 <- acf(sim3.1, lag.max=30,main=expression(''))
pacf1 <- pacf(sim3.1, lag.max=30,main=expression(''))
adf.test(sim3.1)
zroots1 <- polyroot(c(-0.85,1))

dev.off()

#ARIMA(0,0,0)x(1,0,0)12 PHI1=0.85 

sim3.2 <- arima.sim(model = list(ar=c(rep(0,11),-0.85), order=c(12,0,0)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.2,main=expression(''),ylab=' ',xlab=' ')
acf(sim3.2, lag.max=30,main=expression(''))
pacf(sim3.2, lag.max=30,main=expression(''))
adf.test(sim3.2)
zroots2 <- polyroot(c(0.85,rep(0,11),1))
len2 <- abs(zroots2)

#ARIMA(1,0,0)x(0,0,1)12 phi1=-0.8 and THETA1=0.9 

sim3.3 <- arima.sim(model = list(ar=0.8,ma=c(rep(0,11),0.9), order=c(1,0,12)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.3,main=expression(''),ylab=' ',xlab=' ')
acf(sim3.3, lag.max=30,main=expression(''))
pacf(sim3.3, lag.max=30,main=expression(''))
adf.test(sim3.3)
zroots3 <- polyroot(c(0.9,rep(0,11),1))
len3 <- abs(zroots3)

#ARIMA(1,0,0)x(1,0,0)12 phi1=0.7 and PHI1=0.8

sim3.4 <- arima.sim(model = list(ar=c(-0.7,rep(0,10),-0.8, -0.56), order=c(13,0,0)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.4,main=expression(''),ylab=' ',xlab=' ')
acf(sim3.4, lag.max=30,main=expression(''))
pacf(sim3.4, lag.max=30,main=expression(''))
adf.test(sim3.4)
zroots4 <- polyroot(c(0.56, 0.8, rep(0,10), 0.7,1))
len4 <- abs(zroots4)

#ARIMA(2,0,0)x(1,0,0)12 phi1=0.6, phi2=-0.3 and PHI1=0.8

sim3.5 <- arima.sim(model = list(ar=c(-0.6, 0.3,rep(0,9),-0.8, -0.48, 0.24), order=c(14,0,0)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.5,main=expression(''),ylab=' ',xlab=' ')
acf(sim3.5, lag.max=30,main=expression(''))
pacf(sim3.5, lag.max=30,main=expression(''))
adf.test(sim3.5)
zroots5 <- polyroot(c(-0.24,0.48, 0.8, rep(0,9), -0.3,0.6,1))
len5 <- abs(zroots5)

#ARIMA(0,0,1)x(0,0,1)12 theta1=-0.4 and THETA1=0.8 

sim3.6 <- arima.sim(model = list(ma=c(-0.4,rep(0,10),0.8, -0.32), order=c(0,0,13)), n = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(sim3.6,main=expression(''),ylab=' ',xlab=' ')
acf(sim3.6, lag.max=30,main=expression(''))
pacf(sim3.6, lag.max=30,main=expression(''))
adf.test(sim3.6)
