rm(list=ls())

## Predicting Lap Times in a Formula 1 Race
## Using the mathematical formulation for OLS and WLS
## In the end, the relaxation method is used to find the optimal covariance matrix

library(ggplot2)
library(magic)

## load data
load("DataA1.RData")

## Get unique drivers and create mapping
unique_drivers <- sort(unique(AllDat$driverId))
n_drivers <- length(unique_drivers)
driver_map <- setNames(1:n_drivers, unique_drivers)

cat("Number of drivers:", n_drivers, "\n")
cat("Driver IDs:", unique_drivers, "\n")

## Select first driver for analysis
target_driver <- unique_drivers[1]
target_col <- driver_map[as.character(target_driver)]
cat("Target driver:", target_driver, "(column", target_col, ")\n")

## Split data: train = first 51 laps, test = remaining laps
train <- AllDat[AllDat$lap_number <= 51, ]
test <- AllDat[AllDat$lap_number > 51, ]

cat("Train rows:", nrow(train), "\n")
cat("Test rows:", nrow(test), "\n")

## Build design matrix for training data
a_mat <- matrix(0, nrow(train), n_drivers)
for(i in 1:nrow(train)) {
  idx <- driver_map[as.character(train$driverId[i])]
  a_mat[i, idx] <- 1
}

b_mat <- matrix(train$lap_number)
c_mat <- as.integer(train$pitstop)
d_mat <- as.integer(train$pitstop_lagged)
e_mat <- as.integer(train$first_lap)

X <- cbind(a_mat, b_mat, c_mat, d_mat, e_mat)
Y <- train$lap_time

cat("Design matrix dimensions:", dim(X), "\n")

##============================================================================
## OLS ESTIMATION
##============================================================================

theta <- solve(t(X) %*% X, t(X) %*% Y)

n <- length(Y)
p <- dim(X)[2]
err <- Y - X %*% theta
sigma.hat <- sqrt(sum(err^2) / (n - p))

cat("OLS Sigma hat:", sigma.hat, "\n")

## Get parameters for target driver
target_cols <- c(target_col, n_drivers + 1, n_drivers + 2, n_drivers + 3, n_drivers + 4)
theta_target <- theta[target_cols]

## Get training data for target driver
train_target_idx <- train$driverId == target_driver
X_train_target <- X[train_target_idx, target_cols]
Y_train_target <- Y[train_target_idx]

cat("Training rows for target driver:", sum(train_target_idx), "\n")

## Variance of the parameters for target driver
Var.theta_target <- diag(sigma.hat^2 * solve(t(X_train_target) %*% X_train_target))
std.theta_target <- sqrt(Var.theta_target)

## Build design matrix for test data
a_mat2 <- matrix(0, nrow(test), n_drivers)
for(i in 1:nrow(test)) {
  idx <- driver_map[as.character(test$driverId[i])]
  a_mat2[i, idx] <- 1
}

b_mat2 <- matrix(test$lap_number)
c_mat2 <- as.integer(test$pitstop)
d_mat2 <- as.integer(test$pitstop_lagged)
e_mat2 <- as.integer(test$first_lap)

X_test <- cbind(a_mat2, b_mat2, c_mat2, d_mat2, e_mat2)

## Get rows for target driver in test set
test_target_idx <- test$driverId == target_driver
X_test_target <- X_test[test_target_idx, target_cols]

## Predictions for test set
Y_pred_test <- X_test_target %*% theta_target

## Prediction on training set
Y_train_pred <- X_train_target %*% theta_target

## Calculate prediction intervals for test set
Var_pred_test <- diag(sigma.hat^2 * (matrix(1, nrow(X_test_target), nrow(X_test_target)) + 
                                      X_test_target %*% solve(t(X_train_target) %*% X_train_target) %*% t(X_test_target)))

## 95% prediction interval for test
t_val <- qt(0.95, n - p)
interval_pred_test <- cbind(Y_pred_test - t_val * sqrt(Var_pred_test), 
                            Y_pred_test + t_val * sqrt(Var_pred_test))

## 95% prediction interval for training
interval_pred_train <- matrix(0, nrow(X_train_target), 2)
for(i in 1:nrow(X_train_target)) {
  un_root <- 1 + t(X_train_target[i,]) %*% solve(t(X_train_target) %*% X_train_target) %*% X_train_target[i,]
  interval_pred_train[i, 1] <- Y_train_pred[i] - t_val * sigma.hat * sqrt(un_root)
  interval_pred_train[i, 2] <- Y_train_pred[i] + t_val * sigma.hat * sqrt(un_root)
}

## Combine intervals
interval_total_ols <- as.data.frame(rbind(interval_pred_train, interval_pred_test))
Y_fit_ols <- as.data.frame(rbind(Y_train_pred, Y_pred_test))

## Get observed data for target driver
df_target <- AllDat[AllDat$driverId == target_driver, ]

## Plot OLS results
p1 <- ggplot() +
  geom_point(data = df_target, aes(x = lap_number, y = lap_time, color = "Observed \nvalues"), size = 0.8) +
  geom_line(data = Y_fit_ols, aes(x = as.integer(rownames(Y_fit_ols)), y = V1, color = "Fitted \nvalues"), 
            linetype = "solid", size = 0.7) +
  geom_line(data = interval_total_ols, aes(x = as.integer(rownames(interval_total_ols)), y = V1, 
            color = "Prediction \ninterval"), linetype = "twodash", size = 1) +
  geom_line(data = interval_total_ols, aes(x = as.integer(rownames(interval_total_ols)), y = V2, 
            color = "Prediction \ninterval"), linetype = "twodash", size = 1) +
  labs(title = "Fitted and observed values with 95% prediction interval (OLS)", 
       y = "Lap time (ms)", x = "Lap number", color = " ") +
  theme(plot.title = element_text(family = "Helvetica", size = 12, hjust = 0.5), 
        legend.title = element_text(colour = "steelblue", face = "italic", family = "Helvetica"), 
        legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = 8), 
        axis.title = element_text(family = "Helvetica", size = 8, colour = "steelblue4"),
        axis.text = element_text(family = "Courier", colour = "steelblue4", size = 8),
        legend.position = "bottom")

print(p1)

##============================================================================
## WLS - Fixed correlation (rho = 0.5)
##============================================================================

cat("\n--- WLS with fixed correlation ---\n")

## Build covariance matrix with correlation rho = 0.5
rho <- rep(0.5, n_drivers)

## Calculate laps per driver in training set
LapsPerDriver_train <- sapply(unique_drivers, function(d) {
  sum(train$driverId == d)
})

cat("Laps per driver in training:\n")
print(LapsPerDriver_train)

## Initialize Sigma matrix
Sigma <- matrix(0, ncol = 0, nrow = 0)
for(k in seq_along(LapsPerDriver_train)) {
  N <- LapsPerDriver_train[k]
  O <- matrix(data = rep(0:(N-1), N), ncol = N)
  pow <- abs(t(O) - O)
  Block <- rho[k]^pow
  Sigma <- adiag(Sigma, Block)
}

Sigma1 <- Sigma

## Calculate WLS parameters
theta_w <- solve(t(X) %*% solve(Sigma1) %*% X) %*% t(X) %*% solve(Sigma1) %*% Y

## Standard deviation
sigma.hat_w <- sqrt(1/(n-p) * t(Y - X %*% theta_w) %*% solve(Sigma1) %*% (Y - X %*% theta_w))
cat("WLS Sigma hat (fixed rho):", as.numeric(sigma.hat_w), "\n")

## Parameters for target driver
theta_w_target <- theta_w[target_cols]

## Build Sigma matrix for target driver
k <- which(unique_drivers == target_driver)
N <- LapsPerDriver_train[k]
O <- matrix(data = rep(0:(N-1), N), ncol = N)
pow <- abs(t(O) - O)
Sigma_target <- rho[k]^pow

## Variance of parameters
Var.theta_w_target <- diag(as.numeric(sigma.hat_w)^2 * solve(t(X_train_target) %*% solve(Sigma_target) %*% X_train_target))

## Predict lap times
Y_w_target <- X_train_target %*% theta_w_target

## Plot WLS results
data_target_train <- df_target[df_target$lap_number <= 51, ]
Y_w_target_df <- as.data.frame(Y_w_target)

p2 <- ggplot() +
  geom_line(data = data_target_train, aes(x = lap_number, y = lap_time, color = "Observed \nvalues"), 
            linetype = "solid", size = 0.8) +
  geom_line(data = Y_w_target_df, aes(x = as.integer(rownames(Y_w_target_df)), y = V1, 
            color = "Fitted \nvalues"), linetype = "solid", size = 0.8) +
  labs(title = "Fitted and observed values (WLS, fixed rho)", 
       y = "Lap time (ms)", x = "Lap number", color = " ") +
  theme(plot.title = element_text(family = "Helvetica", size = 12, hjust = 0.5), 
        legend.title = element_text(colour = "steelblue", face = "italic", family = "Helvetica"), 
        legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = 8), 
        axis.title = element_text(family = "Helvetica", size = 8, colour = "steelblue4"),
        axis.text = element_text(family = "Courier", colour = "steelblue4", size = 8),
        legend.position = "bottom")

print(p2)

##============================================================================
## WLS - Optimal covariance matrix using relaxation method
##============================================================================

cat("\n--- WLS with relaxation method ---\n")

## Initialize parameters for relaxation
rho_new <- matrix(1, 1, 6)
gamma_first_lap <- matrix(1, 1, 6)
gamma_pit_stop <- matrix(1, 1, 6)
gamma_pit_lag <- matrix(1, 1, 6)

theta_w_opt <- theta_w

## Run relaxation method for 6 iterations
for(i in 1:6) {
  cat("Iteration", i, "\n")
  
  ## Compute residuals
  res <- Y - X %*% theta_w_opt
  
  ## Estimate rho from correlation between residuals at t and t-1
  rho_new[i] <- cor(res[-1], res[-n])
  
  ## Set up Sigma with new rho
  Sigma_3 <- matrix(0, ncol = 0, nrow = 0)
  for(k in seq_along(LapsPerDriver_train)) {
    N <- LapsPerDriver_train[k]
    O <- matrix(data = rep(0:(N-1), N), ncol = N)
    pow <- abs(t(O) - O)
    Block <- rho_new[i]^pow
    Sigma_3 <- adiag(Sigma_3, Block)
  }
  
  ## Calculate standard deviations for different lap types
  sigma_pit <- sd(res[train$pitstop == TRUE])
  sigma_pit_lag <- sd(res[train$pitstop_lagged == 1])
  sigma_first_lap <- sd(res[train$first_lap == 1])
  sigma_denom <- sd(res[train$first_lap != 1 & train$pitstop != TRUE & train$pitstop_lagged != 1])
  
  ## Find gamma factors
  gamma_first_lap[i] <- sigma_first_lap / sigma_denom
  gamma_pit_stop[i] <- sigma_pit / sigma_denom
  gamma_pit_lag[i] <- sigma_pit_lag / sigma_denom
  
  cat("  rho:", rho_new[i], "\n")
  cat("  gamma_first_lap:", gamma_first_lap[i], "\n")
  cat("  gamma_pit_stop:", gamma_pit_stop[i], "\n")
  cat("  gamma_pit_lag:", gamma_pit_lag[i], "\n")
  
  ## Build optimal covariance matrix
  tmp1 <- train$pitstop
  tmp2 <- train$first_lap == 1
  tmp3 <- train$pitstop_lagged == 1
  
  Sigma_4 <- Sigma_3 * 
             (gamma_pit_stop[i]^outer(tmp1, tmp1, FUN = "+")) * 
             (gamma_first_lap[i]^outer(tmp2, tmp2, FUN = "*")) * 
             (gamma_pit_lag[i]^outer(tmp3, tmp3, FUN = "+"))
  
  ## Estimate parameters
  theta_w_opt <- solve(t(X) %*% solve(Sigma_4) %*% X) %*% t(X) %*% solve(Sigma_4) %*% Y
  
  ## Standard deviation
  sigma.hat_w_opt <- as.numeric(sqrt(1/(n-p) * t(Y - X %*% theta_w_opt) %*% solve(Sigma_4) %*% (Y - X %*% theta_w_opt)))
  cat("  sigma hat:", sigma.hat_w_opt, "\n")
}

## Final parameters for target driver
theta_target_final <- theta_w_opt[target_cols]

## Predictions
Y_train_pred_final <- X_train_target %*% theta_target_final
Y_test_pred_final <- X_test_target %*% theta_target_final

## Variance of final parameters
Var.theta_w_final <- diag(sigma.hat_w_opt^2 * solve(t(X) %*% solve(Sigma_4) %*% X))
std.theta_w_final <- sqrt(Var.theta_w_final)

## Extract Sigma for target driver from final Sigma_4
target_train_idx_global <- which(train$driverId == target_driver)
Sigma_target_final <- Sigma_4[target_train_idx_global, target_train_idx_global]

## Combine train and test design matrices for target driver
X_target_combined <- rbind(X_train_target, X_test_target)
Y_target_pred_combined <- rbind(Y_train_pred_final, Y_test_pred_final)

## Calculate prediction intervals
n_total <- nrow(X_target_combined)
interval_pred_final <- matrix(0, n_total, 2)

for(i in 1:n_total) {
  un_root <- 1 + t(X_target_combined[i,]) %*% solve(t(X_train_target) %*% solve(Sigma_target_final) %*% X_train_target) %*% X_target_combined[i,]
  interval_pred_final[i, 1] <- Y_target_pred_combined[i] - t_val * sigma.hat_w_opt * sqrt(un_root)
  interval_pred_final[i, 2] <- Y_target_pred_combined[i] + t_val * sigma.hat_w_opt * sqrt(un_root)
}

interval_pred_final <- as.data.frame(interval_pred_final)
Y_target_pred_combined <- as.data.frame(Y_target_pred_combined)

##============================================================================
## Final plot comparing OLS and WLS
##============================================================================

p_final <- ggplot() + 
  geom_point(data = df_target, aes(x = lap_number, y = lap_time, color = "Observed \nlap times"), size = 0.8) + 
  geom_line(data = interval_pred_final, aes(x = as.integer(rownames(interval_pred_final)), y = V1, 
            color = "Prediction \ninterval"), linetype = "twodash", size = 1) + 
  geom_line(data = interval_pred_final, aes(x = as.integer(rownames(interval_pred_final)), y = V2, 
            color = "Prediction \ninterval"), linetype = "twodash", size = 1) +
  geom_line(data = Y_target_pred_combined, aes(x = as.integer(rownames(Y_target_pred_combined)), y = V1, 
            color = "WLS"), size = 0.7) +
  geom_line(data = Y_fit_ols, aes(x = as.integer(rownames(Y_fit_ols)), y = V1, 
            color = "OLS"), size = 0.6) +
  labs(title = "Fitted values, OLS, observed lap times with 95% prediction interval", 
       y = "Lap time (ms)", x = "Lap number", color = " ") +
  theme(plot.title = element_text(family = "Helvetica", size = 9, hjust = 0.5), 
        legend.title = element_text(colour = "steelblue", face = "italic", family = "Helvetica"), 
        legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = 8), 
        axis.title = element_text(family = "Helvetica", size = 8, colour = "steelblue4"),
        axis.text = element_text(family = "Courier", colour = "steelblue4", size = 8),
        legend.position = "bottom")

print(p_final)

cat("\n=== Analysis complete! ===\n")
cat("Final gamma_first_lap:", gamma_first_lap[6], "\n")
cat("Final gamma_pit_stop:", gamma_pit_stop[6], "\n")
cat("Final gamma_pit_lag:", gamma_pit_lag[6], "\n")
cat("Final rho:", rho_new[6], "\n")