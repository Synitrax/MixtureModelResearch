library(quantmod)
library(mclust)

################################################################################
# GIBBS SAMPLING FUNCTION
################################################################################
run_gibbs <- function(data, x_grid, n_iter = 100, burn_in = 0) {
  n <- length(data)
  y_bar <- mean(data)
  
  # Priors
  mu_0 <- 0; tau_0_sq <- 100 
  a <- 0.01; b <- 0.01        
  
  # Storage for parameters
  mu_samples <- numeric(n_iter)
  sig2_samples <- numeric(n_iter)
  
  # Initial values
  curr_mu <- 0
  curr_sig2 <- 1
  
  # 1. Gibbs Sampling Loop
  for(s in 1:n_iter) {
    # Update Mu (Normal)
    post_prec_mu <- (1/tau_0_sq) + (n/curr_sig2)
    post_mu_mean <- ((mu_0/tau_0_sq) + (n*y_bar/curr_sig2)) / post_prec_mu
    curr_mu <- rnorm(1, post_mu_mean, sqrt(1/post_prec_mu))
    
    # Update Sigma^2 (Inverse-Gamma)
    curr_a <- a + n/2
    curr_b <- b + sum((data - curr_mu)^2)/2
    curr_sig2 <- 1 / rgamma(1, shape = curr_a, rate = curr_b)
    
    mu_samples[s] <- curr_mu
    sig2_samples[s] <- curr_sig2
  }
  
  keep <- (burn_in + 1):n_iter
  mu_final <- mu_samples[keep]
  sig2_final <- sig2_samples[keep]
  S <- length(keep)
  
  phi_matrix <- matrix(NA, nrow = length(x_grid), ncol = S)
  
  for(s in 1:S) {
    phi_matrix[, s] <- dnorm(x_grid, mean = mu_final[s], sd = sqrt(sig2_final[s]))
  }
  
  f_hat <- rowMeans(phi_matrix)
  
  return(list(
    x = x_grid,
    y = f_hat,
    phi_full = phi_matrix,
    mu = mu_final,
    sig2 = sig2_final
  ))
}

################################################################################
# DATA ACQUISITION
################################################################################
tickers = c("AAPL", "UAL", "XOM", "AMZN", "JNJ")
getSymbols(tickers, from = "2022-04-03", to = "2026-04-03")

# Create a list to store CCR data (removing first NA)
stock_list <- list(
  AAPL = na.omit(diff(log(Ad(AAPL)))),
  UAL  = na.omit(diff(log(Ad(UAL)))),
  XOM  = na.omit(diff(log(Ad(XOM)))),
  AMZN = na.omit(diff(log(Ad(AMZN)))),
  JNJ  = na.omit(diff(log(Ad(JNJ))))
)

################################################################################
# COMPARISON LOOP
################################################################################
par(mfrow=c(2,3))

for(name in names(stock_list)) {
  y_real = as.numeric(stock_list[[name]])
  
  # 1. Define Evaluation Grid (xj)
  xj_real = seq(min(y_real), max(y_real), length.out = 100)
  
  # 2. Run Models
  gibbs_out = run_gibbs(y_real, xj_real)
  em_test   = densityMclust(y_real, G=1, verbose = FALSE, plot = FALSE)
  
  # 3. Generate Densities
  # "True" is now the Empirical Density (Kernel Density Estimate)
  true_dens_emp = density(y_real, from=min(xj_real), to=max(xj_real), n=100)$y
  
  bayesian_dens = gibbs_out$y  
  em_dens       = predict(em_test, xj_real)
  
  # 4. Calculate Performance (MSE against Empirical)
  mse_bayesian = mean((bayesian_dens - true_dens_emp)^2)
  mse_em       = mean((em_dens - true_dens_emp)^2)
  
  # 5. Plotting
  plot(density(y_real), main=paste("CCR Analysis:", name), lwd=1, col="gray")
  lines(xj_real, bayesian_dens, col="firebrick", lwd=2)
  lines(xj_real, em_dens, col="royalblue", lwd=2, lty=2)
  
  cat(paste("---", name, "Results ---\n"))
  cat("MSE Bayesian:", mse_bayesian, "\n")
  cat("MSE EM:      ", mse_em, "\n\n")
}

legend("bottomright", legend=c("Empirical", "Bayesian", "EM"), 
       col=c("gray", "firebrick", "royalblue"), lty=c(1,1,2), cex=0.8)
par(mfrow=c(1,1))

################################################################################

# Define the Split Date
split_date <- "2025-04-03"

cat("--- Out-of-Sample Performance (Testing on the Final Year) ---\n\n")

for(name in names(stock_list)) {
  # 1. Get full series and Split Data
  full_series <- stock_list[[name]]
  train_data  <- as.numeric(full_series[index(full_series) < split_date])
  test_data   <- as.numeric(full_series[index(full_series) >= split_date])
  
  # 2. Define Evaluation Grid based on TEST data range
  # This ensures we are evaluating the models where the "future" actually happens
  buffer  <- diff(range(test_data)) * 0.1
  xj_test <- seq(min(test_data) - buffer, max(test_data) + buffer, length.out = 100)
  
  # 3. Train Models
  # Pass the test grid into run_gibbs so it does the math for you
  gibbs_train = run_gibbs(train_data, xj_test)
  em_train    = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  
  # 4. Generate Densities
  # Empirical density of the FUTURE year on the SAME grid
  true_dens_future = density(test_data, from=min(xj_test), to=max(xj_test), n=100)$y
  
  # Pull Bayesian density directly from your function output
  bayesian_dens_pred = gibbs_train$y 
  
  # Predict EM density on the same grid
  em_dens_pred = predict(em_train, xj_test)
  
  # 5. Calculate Out-of-Sample MSE
  mse_bayesian_oos = mean((bayesian_dens_pred - true_dens_future)^2)
  mse_em_oos       = mean((em_dens_pred - true_dens_future)^2)
  
  # 6. Result
  winner <- ifelse(mse_bayesian_oos < mse_em_oos, "Bayesian", "EM")
  
  cat(paste("Stock:", name, "\n"))
  cat("OOS MSE Bayesian:", round(mse_bayesian_oos, 6), "\n")
  cat("OOS MSE EM:      ", round(mse_em_oos, 6), "\n")
  cat("Winner:          ", winner, "\n")
  cat("----------------------------\n")
}
