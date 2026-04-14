library(quantmod)
library(mclust)

################################################################################
# GIBBS SAMPLING FUNCTION
################################################################################
run_gibbs <- function(data, n_iter = 100, burn_in = 0) {
  n = length(data)
  y_bar = mean(data)
  
  mu_0 = 0; tau_0_sq = 100 
  a = 0.01; b = 0.01        
  
  mu_samples = numeric(n_iter)
  sig2_samples = numeric(n_iter)
  
  curr_mu = 0
  curr_sig2 = 1
  
  for(s in 1:n_iter) {
    post_prec_mu = (1/tau_0_sq) + (n/curr_sig2)
    post_mu_mean = ((mu_0/tau_0_sq) + (n*y_bar/curr_sig2)) / post_prec_mu
    curr_mu = rnorm(1, post_mu_mean, sqrt(1/post_prec_mu))
    
    curr_a = a + n/2
    curr_b = b + sum((data - curr_mu)^2)/2
    curr_sig2 = 1 / rgamma(1, shape = curr_a, rate = curr_b)
    
    mu_samples[s] = curr_mu
    sig2_samples[s] = curr_sig2
  }
  
  return(list(mu = mu_samples[(burn_in+1):n_iter], 
              sig2 = sig2_samples[(burn_in+1):n_iter]))
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
par(mfrow=c(2,3)) # Layout for 5 plots + legend space

for(name in names(stock_list)) {
  y_real = as.numeric(stock_list[[name]])
  
  # 1. Define Evaluation Grid (xj)
  xj_real = seq(min(y_real), max(y_real), length.out = 100)
  
  # 2. Run Models
  gibbs_out = run_gibbs(y_real)
  em_test   = densityMclust(y_real, G=1, verbose = FALSE, plot = FALSE)
  
  # 3. Generate Densities
  # "True" is now the Empirical Density (Kernel Density Estimate)
  true_dens_emp = density(y_real, from=min(xj_real), to=max(xj_real), n=1000)$y
  
  bayesian_dens = sapply(xj_real, function(x) mean(dnorm(x, gibbs_out$mu, sqrt(gibbs_out$sig2))))
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










# Define the Split Date
split_date <- "2025-04-03"

cat("--- Out-of-Sample Performance (Testing on the Final Year) ---\n\n")

for(name in names(stock_list)) {
  # Get full series
  full_series <- stock_list[[name]]
  
  # Split Data
  train_data <- as.numeric(full_series[index(full_series) < split_date])
  test_data  <- as.numeric(full_series[index(full_series) >= split_date])
  
  # 1. Train Models on Train Data
  gibbs_train = run_gibbs(train_data)
  em_train    = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  
  # 2. Define Evaluation Grid based on TEST data range
  xj_test = seq(min(test_data), max(test_data), length.out = 1000)
  
  # 3. Generate Densities
  # "True" is now the Empirical Density of the FUTURE year
  true_dens_future = density(test_data, from=min(xj_test), to=max(xj_test), n=1000)$y
  
  # Bayesian: Average of Gaussians using TRAIN samples
  bayesian_dens_pred = sapply(xj_test, function(x) {
    mean(dnorm(x, gibbs_train$mu, sqrt(gibbs_train$sig2)))
  })
  # EM: Single Gaussian using TRAIN parameters
  em_dens_pred = predict(em_train, xj_test)
  
  # 4. Calculate Out-of-Sample MSE
  mse_bayesian_oos = mean((bayesian_dens_pred - true_dens_future)^2)
  mse_em_oos       = mean((em_dens_pred - true_dens_future)^2)
  
  # 5. Result
  winner <- ifelse(mse_bayesian_oos < mse_em_oos, "Bayesian", "EM")
  
  cat(paste("Stock:", name, "\n"))
  cat("OOS MSE Bayesian:", round(mse_bayesian_oos, 4), "\n")
  cat("OOS MSE EM:      ", round(mse_em_oos, 4), "\n")
  cat("Winner:          ", winner, "\n")
  cat("----------------------------\n")
}
