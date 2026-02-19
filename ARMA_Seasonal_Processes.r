# Panagiotis Tsagkaroulis
# Time Series Analysis Script - Runnable Version

rm(list=ls())

# Install and load required packages
if(!require(tseries)) install.packages("tseries", repos='http://cran.rstudio.com/')
if(!require(tidyr)) install.packages("tidyr", repos='http://cran.rstudio.com/')
if(!require(ggplot2)) install.packages("ggplot2", repos='http://cran.rstudio.com/')

library(tseries)
library(tidyr)
library(ggplot2)

##============================================================================
## Part 1: Check stationarity of MA(3) process
##============================================================================

cat("\n=== Part 1: MA(3) Process Analysis ===\n")

# Check stationarity by finding roots
zroots <- polyroot(c(1, 1, 1, 1))
cat("Roots of MA polynomial:\n")
print(zroots)
cat("Absolute values of roots:\n")
print(abs(zroots))

##============================================================================
## Part 2: Simulate 10 realizations with 200 observations
##============================================================================

cat("\n=== Part 2: Simulating MA(3) Process ===\n")

set.seed(123)  # For reproducibility

# Simulate 10 realizations
sim1 <- cbind("Time" = 1:200, 
              replicate(10, arima.sim(model = list(ma = c(1, 1, 1), order = c(0, 0, 3)), 
                                     sd = 0.1, n = 200)))
colnames(sim1) <- c('Time', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10')

# Convert to long format for plotting
dfsim1 <- pivot_longer(as.data.frame(sim1), cols = 2:11, names_to = "name", values_to = "value")

# Plot the 10 realizations
p <- ggplot(dfsim1, aes(y = value, color = name, x = Time)) + 
  geom_line() +
  theme_light() +
  theme(axis.title = element_text(size = 10), 
        legend.position = "none",
        axis.title.x = element_blank()) +
  labs(title = "10 Realizations of MA(3) Process", y = "Value")

print(p)

##============================================================================
## Part 3: Estimate and plot ACF for each realization
##============================================================================

cat("\n=== Part 3: ACF Analysis ===\n")

# Calculate ACF for each realization
acf_list <- lapply(2:11, function(i) {
  acf(sim1[, i], lag.max = 24, plot = FALSE)$acf
})

acfsim1 <- cbind("Lag" = 0:24, do.call(cbind, acf_list))
colnames(acfsim1) <- c("Lag", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

dfacfsim1 <- pivot_longer(as.data.frame(acfsim1), cols = 2:11, names_to = "name", values_to = "value")

# Plot ACF as lines
p1 <- ggplot(dfacfsim1, aes(y = value, color = name, x = Lag)) + 
  geom_line() +
  theme_light() +
  theme(axis.title = element_text(size = 10), 
        legend.position = "none",
        axis.title.x = element_blank()) + 
  labs(y = "ACF", title = "ACF for 10 Realizations") +
  geom_hline(yintercept = 0, linetype = "dashed")

print(p1)

# Alternative: Plot ACF as segments (classical form)
p1_seg <- ggplot(dfacfsim1, aes(y = value, color = name, x = Lag)) + 
  geom_segment(aes(xend = Lag, yend = 0)) +
  theme_light() +
  theme(axis.title = element_text(size = 10), 
        legend.position = "none") + 
  labs(y = "ACF", x = "Lag", title = "ACF for 10 Realizations (Segment Plot)") +
  geom_hline(yintercept = 0)

print(p1_seg)

##============================================================================
## Part 4: Estimate and plot PACF for each realization
##============================================================================

cat("\n=== Part 4: PACF Analysis ===\n")

# Calculate PACF for each realization
pacf_list <- lapply(2:11, function(i) {
  pacf(sim1[, i], lag.max = 24, plot = FALSE)$acf
})

pacfsim1 <- cbind("Lag" = 1:24, do.call(cbind, pacf_list))
colnames(pacfsim1) <- c("Lag", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

dfpacfsim1 <- pivot_longer(as.data.frame(pacfsim1), cols = 2:11, names_to = "name", values_to = "value")

# Plot PACF as lines
p2 <- ggplot(dfpacfsim1, aes(y = value, color = name, x = Lag)) + 
  geom_line() +
  theme_light() +
  theme(axis.title = element_text(size = 10), 
        legend.position = "none",
        axis.title.x = element_blank()) + 
  labs(y = "PACF", title = "PACF for 10 Realizations") +
  geom_hline(yintercept = 0, linetype = "dashed")

print(p2)

# Alternative: Plot PACF as segments
p2_seg <- ggplot(dfpacfsim1, aes(y = value, color = name, x = Lag)) + 
  geom_segment(aes(xend = Lag, yend = 0)) +
  theme_light() +
  theme(axis.title = element_text(size = 10), 
        legend.position = "none") + 
  labs(y = "PACF", x = "Lag", title = "PACF for 10 Realizations (Segment Plot)") +
  geom_hline(yintercept = 0)

print(p2_seg)

##============================================================================
## Part 5: Calculate variance of each realization
##============================================================================

cat("\n=== Part 5: Variance Analysis ===\n")

Variance <- apply(sim1[, 2:11], 2, var)
cat("Variance of each realization:\n")
print(Variance)
cat("\nMean variance:", mean(Variance), "\n")

##============================================================================
## Part 6: Predicting temperatures in district heating network
##============================================================================

cat("\n=== Part 6: Temperature Prediction ===\n")

# Create data frame
year <- c(rep(2016, 4), rep(2017, 11))
month <- c(9:12, 1:11)
temp <- c(46.6, 49.5, 60.3, 59.2, 59.5, 61.9, 59.7, 60.1, 57.8, 49.7, 49.7, 50.1, 48.6, 54.5, 62.3)
df <- data.frame(year, month, temp)

cat("Temperature data:\n")
print(df)

# Transform data (mu = 55)
X.t <- df$temp - 55
X.len <- length(X.t)

# Predictions
X.hat2017M12 <- 0.5 * X.t[X.len] - 0.3 * X.t[X.len - 1] + 0.9 * X.t[X.len - 11] - 
                0.45 * X.t[X.len - 12] + 0.27 * X.t[X.len - 13]

X.hat2018M1 <- 0.5 * X.hat2017M12 - 0.3 * X.t[X.len] + 0.9 * X.t[X.len - 10] - 
               0.45 * X.t[X.len - 11] + 0.27 * X.t[X.len - 12]

# Back-transformation
Y.hat2017M12 <- X.hat2017M12 + 55
Y.hat2018M1 <- X.hat2018M1 + 55

cat("\nPredicted temperature for 2017 M12:", Y.hat2017M12, "\n")
cat("Predicted temperature for 2018 M1:", Y.hat2018M1, "\n")

# Variance of prediction error
sigma.epsilon <- 0.5
V.epsilon2017M12 <- sigma.epsilon^2
V.epsilon2018M1 <- (1 + (0.5)^2) * (sigma.epsilon^2)

# 95% confidence interval
inter <- qnorm(0.975)

Y2017M12.low <- Y.hat2017M12 - inter * sqrt(V.epsilon2017M12)
Y2017M12.high <- Y.hat2017M12 + inter * sqrt(V.epsilon2017M12)
Y2018M1.low <- Y.hat2018M1 - inter * sqrt(V.epsilon2018M1)
Y2018M1.high <- Y.hat2018M1 + inter * sqrt(V.epsilon2018M1)

cat("\n95% CI for 2017 M12: [", Y2017M12.low, ",", Y2017M12.high, "]\n")
cat("95% CI for 2018 M1: [", Y2018M1.low, ",", Y2018M1.high, "]\n")

# Plot predictions with confidence intervals
prediction.interval.low <- c(Y2017M12.low, Y2018M1.low)
prediction.interval.high <- c(Y2017M12.high, Y2018M1.high)

time <- 1:15
plot(time, df$temp, 
     xlim = c(0, 18), 
     ylim = c(min(df$temp) - 5, max(df$temp) + 5),
     main = '95% Confidence Interval for Temperature Predictions',
     xlab = 'Time', 
     ylab = 'Temperature',
     pch = 16)

points(16, Y.hat2017M12, col = "red", pch = 16) 
points(17, Y.hat2018M1, col = "red", pch = 16)
points(16:17, prediction.interval.low, col = 2, pch = 95)
points(16:17, prediction.interval.high, col = 2, pch = 95)
lines(time, df$temp)

par(xpd = TRUE)
legend(x = "bottomright", 
       inset = 0.01, 
       legend = c("Observed", "Prediction", "95% CI"), 
       pch = c(16, 16, 95), 
       col = c(1, 2, 2), 
       cex = 0.8)

##============================================================================
## Part 7: Simulating Seasonal ARIMA Processes
##============================================================================

cat("\n=== Part 7: Seasonal ARIMA Simulations ===\n")

# Custom ARIMA simulation function (extended version)
arima.sim2 <- function(model, n, rand.gen = rnorm, innov = rand.gen(n, ...), 
                       n.start = NA, start.innov = rand.gen(n.start, ...), ...) {
  if (!is.list(model)) 
    stop("'model' must be list")
  if (n <= 0L) 
    stop("'n' must be strictly positive")
  
  p <- length(model$ar)
  if (p) {
    minroots <- min(Mod(polyroot(c(1, -model$ar))))
    if (minroots <= 1) 
      warning("'ar' part of model is not stationary")
  }
  
  q <- length(model$ma)
  if (is.na(n.start)) 
    n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)), 0)
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
  
  x <- ts(c(start.innov[seq_len(n.start)], innov[1L:n]), start = 1 - n.start)
  
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

##----------------------------------------------------------------------------
## Model 1: ARIMA(1,0,0)x(0,0,0)12 with phi1 = 0.85
##----------------------------------------------------------------------------

cat("\n--- Model 1: AR(1) with phi1 = 0.85 ---\n")

set.seed(123)
sim3.1 <- arima.sim(model = list(ar = 0.85, order = c(1, 0, 0)), n = 200)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.1, main = 'ARIMA(1,0,0) - AR(1) Process', ylab = 'Value', xlab = 'Time')
acf(sim3.1, lag.max = 30, main = 'ACF')
pacf(sim3.1, lag.max = 30, main = 'PACF')

adf_result <- adf.test(sim3.1)
cat("ADF test p-value:", adf_result$p.value, "\n")

zroots1 <- polyroot(c(1, -0.85))
cat("Roots:", zroots1, "\n")
cat("Absolute values:", abs(zroots1), "\n")

dev.off()

##----------------------------------------------------------------------------
## Model 2: ARIMA(0,0,0)x(1,0,0)12 with PHI1 = -0.85
##----------------------------------------------------------------------------

cat("\n--- Model 2: Seasonal AR(1) with PHI1 = -0.85 ---\n")

set.seed(124)
sim3.2 <- arima.sim(model = list(ar = c(rep(0, 11), -0.85), order = c(12, 0, 0)), n = 200)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.2, main = 'ARIMA(0,0,0)x(1,0,0)12 - Seasonal AR', ylab = 'Value', xlab = 'Time')
acf(sim3.2, lag.max = 30, main = 'ACF')
pacf(sim3.2, lag.max = 30, main = 'PACF')

adf_result2 <- adf.test(sim3.2)
cat("ADF test p-value:", adf_result2$p.value, "\n")

zroots2 <- polyroot(c(1, rep(0, 11), -0.85))
len2 <- abs(zroots2)
cat("Absolute values of roots:", len2, "\n")

dev.off()

##----------------------------------------------------------------------------
## Model 3: ARIMA(1,0,0)x(0,0,1)12 with phi1 = 0.8, THETA1 = 0.9
##----------------------------------------------------------------------------

cat("\n--- Model 3: ARIMA(1,0,0)x(0,0,1)12 ---\n")

set.seed(125)
sim3.3 <- arima.sim(model = list(ar = 0.8, ma = c(rep(0, 11), 0.9), order = c(1, 0, 12)), n = 200)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.3, main = 'ARIMA(1,0,0)x(0,0,1)12', ylab = 'Value', xlab = 'Time')
acf(sim3.3, lag.max = 30, main = 'ACF')
pacf(sim3.3, lag.max = 30, main = 'PACF')

adf_result3 <- adf.test(sim3.3)
cat("ADF test p-value:", adf_result3$p.value, "\n")

zroots3 <- polyroot(c(1, rep(0, 11), 0.9))
len3 <- abs(zroots3)
cat("Absolute values of MA roots:", len3, "\n")

dev.off()

##----------------------------------------------------------------------------
## Model 4: ARIMA(1,0,0)x(1,0,0)12 with phi1 = 0.7, PHI1 = 0.8
##----------------------------------------------------------------------------

cat("\n--- Model 4: ARIMA(1,0,0)x(1,0,0)12 ---\n")

set.seed(126)
# For ARIMA(1,0,0)x(1,0,0)12: (1-phi1*B)(1-PHI1*B^12) = 1 - phi1*B - PHI1*B^12 + phi1*PHI1*B^13
# With phi1=0.7, PHI1=0.8: coefficients are [0.7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8, -0.56]
tryCatch({
  sim3.4 <- arima.sim(model = list(ar = c(0.7, rep(0, 10), 0.8, -0.56), order = c(13, 0, 0)), n = 200)
}, warning = function(w) {
  cat("Warning in model simulation:", conditionMessage(w), "\n")
  # Continue anyway
  sim3.4 <<- arima.sim(model = list(ar = c(0.7, rep(0, 10), 0.8, -0.56), order = c(13, 0, 0)), n = 200)
}, error = function(e) {
  cat("Error in model simulation:", conditionMessage(e), "\n")
  cat("Using a more conservative AR specification...\n")
  # Use more conservative parameters
  sim3.4 <<- arima.sim(model = list(ar = c(0.5, rep(0, 10), 0.5, -0.25), order = c(13, 0, 0)), n = 200)
})

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.4, main = 'ARIMA(1,0,0)x(1,0,0)12', ylab = 'Value', xlab = 'Time')
acf(sim3.4, lag.max = 30, main = 'ACF')
pacf(sim3.4, lag.max = 30, main = 'PACF')

adf_result4 <- adf.test(sim3.4)
cat("ADF test p-value:", adf_result4$p.value, "\n")

# Check AR polynomial roots
zroots4 <- polyroot(c(1, -0.7, rep(0, 10), -0.8, 0.56))
len4 <- abs(zroots4)
cat("Absolute values of AR roots (should be > 1 for stationarity):\n")
print(len4)

dev.off()

##----------------------------------------------------------------------------
## Model 5: ARIMA(2,0,0)x(1,0,0)12 with phi1 = 0.6, phi2 = -0.3, PHI1 = 0.8
##----------------------------------------------------------------------------

cat("\n--- Model 5: ARIMA(2,0,0)x(1,0,0)12 ---\n")

set.seed(127)
# For ARIMA(2,0,0)x(1,0,0)12: (1-phi1*B-phi2*B^2)(1-PHI1*B^12)
# With phi1=0.6, phi2=-0.3, PHI1=0.8
# Coefficients: [0.6, -0.3, 0, ..., 0, 0.8, -0.48, 0.24]
tryCatch({
  sim3.5 <- arima.sim(model = list(ar = c(0.6, -0.3, rep(0, 9), 0.8, -0.48, 0.24), 
                                   order = c(14, 0, 0)), n = 200)
}, warning = function(w) {
  cat("Warning in model simulation:", conditionMessage(w), "\n")
  # Continue anyway
  sim3.5 <<- arima.sim(model = list(ar = c(0.6, -0.3, rep(0, 9), 0.8, -0.48, 0.24), 
                                    order = c(14, 0, 0)), n = 200)
}, error = function(e) {
  cat("Error in model simulation:", conditionMessage(e), "\n")
  cat("Using a more conservative AR specification...\n")
  # Use more conservative parameters
  sim3.5 <<- arima.sim(model = list(ar = c(0.5, -0.2, rep(0, 9), 0.6, -0.3, 0.12), 
                                    order = c(14, 0, 0)), n = 200)
})

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.5, main = 'ARIMA(2,0,0)x(1,0,0)12', ylab = 'Value', xlab = 'Time')
acf(sim3.5, lag.max = 30, main = 'ACF')
pacf(sim3.5, lag.max = 30, main = 'PACF')

adf_result5 <- adf.test(sim3.5)
cat("ADF test p-value:", adf_result5$p.value, "\n")

# Check AR polynomial roots
zroots5 <- polyroot(c(1, -0.6, 0.3, rep(0, 9), -0.8, 0.48, -0.24))
len5 <- abs(zroots5)
cat("Absolute values of AR roots (should be > 1 for stationarity):\n")
print(len5)

dev.off()

##----------------------------------------------------------------------------
## Model 6: ARIMA(0,0,1)x(0,0,1)12 with theta1 = -0.4, THETA1 = 0.8
##----------------------------------------------------------------------------

cat("\n--- Model 6: ARIMA(0,0,1)x(0,0,1)12 ---\n")

set.seed(128)
sim3.6 <- arima.sim(model = list(ma = c(-0.4, rep(0, 10), 0.8, -0.32), order = c(0, 0, 13)), n = 200)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
plot(sim3.6, main = 'ARIMA(0,0,1)x(0,0,1)12', ylab = 'Value', xlab = 'Time')
acf(sim3.6, lag.max = 30, main = 'ACF')
pacf(sim3.6, lag.max = 30, main = 'PACF')

adf_result6 <- adf.test(sim3.6)
cat("ADF test p-value:", adf_result6$p.value, "\n")

# Check MA polynomial roots (note: for invertibility, roots should be > 1)
zroots6 <- polyroot(c(1, -0.4, rep(0, 10), 0.8, -0.32))
len6 <- abs(zroots6)
cat("Absolute values of MA roots (should be > 1 for invertibility):\n")
print(len6)

dev.off()

cat("\n=== Script completed successfully! ===\n")