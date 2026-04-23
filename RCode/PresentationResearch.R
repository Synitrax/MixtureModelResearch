# -------------------------------------------------------------------------
# FUNCTION: NVMunmix
# Theory & Algorithm: Dr. Hasan Hamdan
# R Implementation & Optimization: Nathan Carter
# -------------------------------------------------------------------------
source("c:/qpsolve/qpsolve.R")


install.packages("mclust")
install.packages("quadprog")
install.packages("zoo")


library(mclust)
library(zoo)
library(quadprog)


NVMunmix = function(x,b=.4,rgrid,mugrid,tol){
a = -max(abs(x))-(b*2)
cc = max(abs(x))+(b*2)
r= length(rgrid) 
m=length(mugrid) 
#denx = density(x,width=b,from=a,to=cc)
denx = density(x)
# this estimates the density of x, larger the n, better the estimate, slower the computation
xgrid = denx$x # initialize x grid, using the same x-points used by the density function
q=length(xgrid)

phi = matrix(0.0, nrow=q ,ncol=r*m) # nrow = the n in the denx line.
for (h in 1:m) # h counts the u's
  {for (j in 1:r) # j counts the r's
for (i in 1:q) # why 11 again? # i counts the x's
     phi[i,(h-1)*r + j] = dnorm(xgrid[i],mugrid[h],rgrid[j]) } 


 y = denx$y # height of each choosen x value
p1=qpsolve( y, phi, pi.total=1.0 )$pi # feeds super vector to qpsolve program, they have the wieghts for each combination
p1 = c(p1) # combines p1 into vector
 xxx=matrix(nrow = r*m, ncol = 3)  

 for (i in 1:m)
{
   for (j in 1:r)
    xxx[j+(i-1)*r, ] = c(mugrid[i], rgrid[j], p1[j+(i-1)*r]) 
}



xxx=data.frame(xxx)
xxx=split(xxx,xxx[,3]>tol)$'TRUE'
w=sum(xxx[,3])
xxx[,3]=xxx[,3]/w # This is to normalize the remaining terms

Ymix=c(1:q)
for(i in 1:q)

Ymix[i]=sum(xxx[,3]*dnorm(xgrid[i],xxx[,1],xxx[,2]))
x=denx$x
#YBayes =.5966*dnorm(x,536.7,1/sqrt(.07203))+.4034*dnorm(x,548.8,1/sqrt(.07203))
YBayes =0
dd=max(c(Ymix, y,YBayes))
#plot(denx$x,denx$y,type="l",col="red", ylim=c(0,dd))
#lines(denx$x, Ymix,col="blue", lty=4)

#lines(x,YBayes, lty=8)

list(xxx, w, dim(xxx))

}

run_gibbs <- function(data, x_grid, M = 2, n_iter = 100) {
  n <- length(data)
  
  # Storage: Dimensions are [Iteration, Component]
  mu_samples <- matrix(0, nrow = n_iter, ncol = M)
  sig2_samples <- matrix(1, nrow = n_iter, ncol = M)
  pi_samples <- matrix(1/M, nrow = n_iter, ncol = M) # Weights
  
  # Initial values
  curr_mu <- seq(min(data), max(data), length.out = M)
  curr_sig2 <- rep(var(data), M)
  curr_pi <- rep(1/M, M)
  
  for(s in 1:n_iter) {
    # STEP 1: Assign data to components (Latent Variables Z)
    # This calculates the probability each xi belongs to component m
    probs <- matrix(0, n, M)
    for(m in 1:M) {
      probs[, m] <- curr_pi[m] * dnorm(data, curr_mu[m], sqrt(curr_sig2[m]))
    }
    probs <- probs / rowSums(probs)
    # Randomly assign each point to a component based on probs
    z <- apply(probs, 1, function(p) sample(1:M, size = 1, prob = p))
    
    # STEP 2: Update parameters for each component m
    for(m in 1:M) {
      data_m <- data[z == m]
      n_m <- length(data_m)
      
      if(n_m > 0) {
        # Update Mu_m (Normal)
        y_bar_m <- mean(data_m)
        post_prec <- (1/100) + (n_m/curr_sig2[m])
        post_mean <- ((0/100) + (n_m*y_bar_m/curr_sig2[m])) / post_prec
        curr_mu[m] <- rnorm(1, post_mean, sqrt(1/post_prec))
        
        # Update Sigma2_m (Inverse-Gamma)
        curr_a <- 0.01 + n_m/2
        curr_b <- 0.01 + sum((data_m - curr_mu[m])^2)/2
        curr_sig2[m] <- 1 / rgamma(1, shape = curr_a, rate = curr_b)
      }
    }
    
    # STEP 3: Update weights (Dirichlet -> simplified for 2 components)  
    mu_samples[s, ] <- curr_mu
    sig2_samples[s, ] <- curr_sig2
  }
  
  # STEP 4: Density Estimation
  # f_hat(x) = 1/S * Sum_s [ Sum_m (pi_m * phi(x | mu_sm, sig_sm)) ]
  f_grid_samples <- matrix(0, nrow = length(x_grid), ncol = n_iter)
  for(s in 1:n_iter) {
    comp_densities <- matrix(0, length(x_grid), M)
    for(m in 1:M) {
      comp_densities[,m] <- (1/M) * dnorm(x_grid, mu_samples[s,m], sqrt(sig2_samples[s,m]))
    }
    f_grid_samples[, s] <- rowSums(comp_densities)
  }
  
  f_hat <- rowMeans(f_grid_samples)
  return(list(x = x_grid, y = f_hat))
}

# Read file, explicitly telling R NOT to turn strings into factors
stock_data_raw <- read.csv("stock_data.csv", row.names = 1, stringsAsFactors = FALSE)

# If the first column is named "Index", remove it so only prices remain
if("Index" %in% colnames(stock_data_raw)) {
  stock_data_raw$Index <- NULL
}

# Ensure all columns are numeric (this fixes the 'factor' error)
stock_data_raw[] <- lapply(stock_data_raw, function(x) as.numeric(as.character(x)))

# Convert to a list of CCR (Log Returns)
stock_list <- lapply(stock_data_raw, function(x) {
  prices <- na.omit(x)
  # diff(log(x)) calculates: log(P_t) - log(P_t-1)
  returns <- diff(log(prices))
  return(returns)
})

# Get the dates from the row names
# We remove the first date because diff() results in n-1 observations
all_dates <- as.Date(rownames(stock_data_raw))
returns_dates <- all_dates[-1] 

# Quick check:
if(exists("stock_list")) {
  cat("Data loaded successfully. Stocks identified:", names(stock_list), "\n")
}

#########################################################################################################
dev.new() 
par(mfrow = c(2, 3))

for(name in names(stock_list)) {
  # --- 1. DATA PREP ---
  train_data <- stock_list[[name]]
  test_data  <- train_data 
  xj_test = seq(min(test_data), max(test_data), length.out = 1000)
  
  # --- 2. MODEL FITS ---
  # Bayesian (Gibbs)
  gibbs_train = run_gibbs(train_data, xj_test)
  bayesian_pred = gibbs_train$y
  
  # EM (Mclust)
  em_train = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  em_pred  = predict(em_train, xj_test)
  
  # NVMunmix
  mu_grid <- seq(min(train_data), max(train_data), length.out = 25)
  r_grid  <- seq(sd(train_data)*0.3, sd(train_data)*1.2, length.out = 15)
  nvm_res <- tryCatch({
    NVMunmix(train_data, b=0.4, rgrid=r_grid, mugrid=mu_grid, tol=0.001)
  }, error = function(e) return(NULL))
  
  if(!is.null(nvm_res)) {
    mix_df <- nvm_res[[1]]
    nvm_pred <- sapply(xj_test, function(x) sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2])))
  } else {
    nvm_pred <- rep(0, length(xj_test))
  }

  # --- 3. PLOTTING ---
  plot(density(test_data), main=paste("In-Sample Fit:", name), lwd=1, col="gray")
  lines(xj_test, bayesian_pred, col="firebrick", lwd=2)
  lines(xj_test, em_pred, col="royalblue", lwd=2, lty=2)
  lines(xj_test, nvm_pred, col="darkgreen", lwd=2, lty=3)
  legend("topright", legend=c("Empirical", "Bayesian", "EM", "NVM"), 
         col=c("gray", "firebrick", "royalblue", "darkgreen"), lty=c(1,1,2,3), cex=0.7)

  # --- 4. EVALUATION: MSE ---
  true_dens_in = density(test_data, from=min(xj_test), to=max(xj_test), n=1000)$y
  mse_bayesian = mean((bayesian_pred - true_dens_in)^2)
  mse_em       = mean((em_pred - true_dens_in)^2)
  mse_nvm      = mean((nvm_pred - true_dens_in)^2)

  # --- 5. EVALUATION: LOG-LIKELIHOOD ---
  ll_bayesian <- sum(log(approx(xj_test, bayesian_pred, xout = test_data, rule=2)$y + 1e-10))
  ll_em       <- sum(log(predict(em_train, test_data) + 1e-10))
  
  nvm_dens_at_data <- sapply(test_data, function(x) {
    if(!is.null(nvm_res)) sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2])) else 1e-10
  })
  ll_nvm <- sum(log(nvm_dens_at_data + 1e-10))

  # --- 6. EVALUATION: K-S STATISTIC ---
  ecd_func   <- ecdf(test_data)
  actual_cdf <- ecd_func(xj_test)
  dx         <- xj_test[2] - xj_test[1]
  
  cdf_bayesian <- cumsum(bayesian_pred * dx); cdf_bayesian <- cdf_bayesian / (max(cdf_bayesian) + 1e-10)
  cdf_em       <- cumsum(em_pred * dx);       cdf_em       <- cdf_em / (max(cdf_em) + 1e-10)
  cdf_nvm      <- cumsum(nvm_pred * dx);      cdf_nvm      <- cdf_nvm / (max(cdf_nvm) + 1e-10)
  
  ks_bayesian <- max(abs(actual_cdf - cdf_bayesian))
  ks_em       <- max(abs(actual_cdf - cdf_em))
  ks_nvm      <- max(abs(actual_cdf - cdf_nvm))

  # --- 8. EVALUATION: CHI-SQUARE ---
  n_obs <- length(test_data)
  k <- floor(sqrt(n_obs)) # Number of bins
  breaks <- quantile(test_data, probs = seq(0, 1, length.out = k + 1))
  observed_counts <- as.vector(table(cut(test_data, breaks = breaks, include.lowest = TRUE)))

  get_expected <- function(pred_y, grid_x, breaks, n_total) {
    probs <- sapply(1:(length(breaks)-1), function(i) {
      bin_mask <- grid_x >= breaks[i] & grid_x <= breaks[i+1]
      if(sum(bin_mask) < 2) return(1e-10)
      bin_x <- grid_x[bin_mask]; bin_y <- pred_y[bin_mask]
      sum(diff(bin_x) * (bin_y[-1] + bin_y[-length(bin_y)]) / 2)
    })
    probs <- probs / (sum(probs) + 1e-10)
    return(probs * n_total)
  }

  # --- 9. EVALUATION: BIC ---
  n_obs <- length(test_data)
  
  # Define degrees of freedom (k) for each model
  k_bayesian <- 2  # Mean and Variance
  k_em       <- 2  # G=1: Mean and Variance
  # For NVM, k is the number of components with significant weight
  k_nvm      <- if(!is.null(nvm_res)) sum(mix_df[,3] > 1e-4) * 2 else 1 

  bic_bayesian <- k_bayesian * log(n_obs) - 2 * ll_bayesian
  bic_em       <- k_em * log(n_obs) - 2 * ll_em
  bic_nvm      <- k_nvm * log(n_obs) - 2 * ll_nvm
  exp_bayesian <- get_expected(bayesian_pred, xj_test, breaks, n_obs)
  exp_em       <- get_expected(em_pred, xj_test, breaks, n_obs)
  exp_nvm      <- get_expected(nvm_pred, xj_test, breaks, n_obs)

  chisq_bayesian <- sum((observed_counts - exp_bayesian)^2 / (exp_bayesian + 1e-10))
  chisq_em       <- sum((observed_counts - exp_em)^2 / (exp_em + 1e-10))
  chisq_nvm      <- sum((observed_counts - exp_nvm)^2 / (exp_nvm + 1e-10))

  # --- 7. PRINT SUMMARY TABLE ---
  cat(paste("\nStock Analysis:", name, "(IN-SAMPLE)\n"))
  cat(sprintf("%-15s | %-12s | %-12s | %-12s\n", "Metric", "Bayesian", "EM", "NVM"))
  cat(rep("-", 60), "\n")
  cat(sprintf("%-15s | %-12.6f | %-12.6f | %-12.6f\n", "MSE (Lower)", mse_bayesian, mse_em, mse_nvm))
  cat(sprintf("%-15s | %-12.2f | %-12.2f | %-12.2f\n", "Log-Lik (Higher)", ll_bayesian, ll_em, ll_nvm))
  cat(sprintf("%-15s | %-12.4f | %-12.4f | %-12.4f\n", "K-S (Lower)", ks_bayesian, ks_em, ks_nvm))
  cat(sprintf("%-15s | %-12.4f | %-12.4f | %-12.4f\n", "Chi-Sq (Lower)", chisq_bayesian, chisq_em, chisq_nvm))
  cat(sprintf("%-15s | %-12.2f | %-12.2f | %-12.2f\n", "BIC (Lower)", bic_bayesian, bic_em, bic_nvm))
  
  mses <- c(Bayesian = mse_bayesian, EM = mse_em, NVM = mse_nvm)
  chisqs <- c(Bayesian = chisq_bayesian, EM = chisq_em, NVM = chisq_nvm)
  bics <- c(Bayesian = bic_bayesian, EM = bic_em, NVM = bic_nvm)
  cat("Winner (MSE):   ", names(which.min(mses)), "\n")
  cat("Winner (Chi-Sq):", names(which.min(chisqs)), "\n")
  cat("Winner (BIC):   ", names(which.min(bics)), "\n")
  cat("------------------------------------------------------------\n")
}
############################################################################
dev.new() 
par(mfrow = c(2, 3))

split_date <- as.Date("2025-04-03")

cat("--- Out-of-Sample Performance: Bayesian vs EM vs NVM ---\n\n")

for(name in names(stock_list)) {
  # --- 1. DATA PREP & SPLIT ---
  full_series <- stock_list[[name]]
  
  # Split based on the dates vector
  train_data <- full_series[returns_dates < split_date]
  test_data  <- full_series[returns_dates >= split_date]
  
  # Define Evaluation Grid based on TEST data range (where we evaluate performance)
  xj_test = seq(min(test_data), max(test_data), length.out = 1000)
  
  # --- 2. MODEL TRAINING (on train_data) & PREDICTION (on xj_test) ---
  
  # Bayesian (Gibbs)
  gibbs_train = run_gibbs(train_data, xj_test)
  bayesian_pred = gibbs_train$y
  
  # EM (Mclust) - Fit on train, predict on test grid
  em_train = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  em_pred  = predict(em_train, xj_test)
  
  # NVMunmix
  mu_grid <- seq(min(train_data), max(train_data), length.out = 20)
  r_grid  <- seq(sd(train_data)*0.5, sd(train_data)*1.5, length.out = 10)
  
  nvm_res <- tryCatch({
    NVMunmix(train_data, b=0.4, rgrid=r_grid, mugrid=mu_grid, tol=0.001)
  }, error = function(e) return(NULL))
  
  if(!is.null(nvm_res)) {
    mix_df <- nvm_res[[1]]
    nvm_pred <- sapply(xj_test, function(x) sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2])))
  } else {
    nvm_pred <- rep(0, length(xj_test))
  }

  # --- 3. PLOTTING ---
  plot(density(test_data), main=paste("OOS Predict:", name), lwd=1, col="gray")
  lines(xj_test, bayesian_pred, col="firebrick", lwd=2)
  lines(xj_test, em_pred, col="royalblue", lwd=2, lty=2)
  lines(xj_test, nvm_pred, col="darkgreen", lwd=2, lty=3)
  legend("topright", legend=c("True Test Density", "Bayesian", "EM", "NVM"), 
         col=c("gray", "firebrick", "royalblue", "darkgreen"), lty=c(1,1,2,3), cex=0.6)

  # --- 4. EVALUATION: MSE ---
  true_dens_future = density(test_data, from=min(xj_test), to=max(xj_test), n=1000)$y
  mse_bayesian = mean((bayesian_pred - true_dens_future)^2)
  mse_em       = mean((em_pred - true_dens_future)^2)
  mse_nvm      = mean((nvm_pred - true_dens_future)^2)

  # --- 5. EVALUATION: LOG-LIKELIHOOD (on test_data) ---
  ll_bayesian <- sum(log(approx(xj_test, bayesian_pred, xout = test_data, rule=2)$y + 1e-10))
  ll_em       <- sum(log(predict(em_train, test_data) + 1e-10))
  
  nvm_dens_at_data <- sapply(test_data, function(x) {
    if(!is.null(nvm_res)) sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2])) else 1e-10
  })
  ll_nvm <- sum(log(nvm_dens_at_data + 1e-10))

  # --- 6. EVALUATION: K-S STATISTIC ---
  ecd_func   <- ecdf(test_data)
  actual_cdf <- ecd_func(xj_test)
  dx         <- xj_test[2] - xj_test[1]
  
  cdf_bayesian <- cumsum(bayesian_pred * dx); cdf_bayesian <- cdf_bayesian / (max(cdf_bayesian) + 1e-10)
  cdf_em       <- cumsum(em_pred * dx);       cdf_em       <- cdf_em / (max(cdf_em) + 1e-10)
  cdf_nvm      <- cumsum(nvm_pred * dx);      cdf_nvm      <- cdf_nvm / (max(cdf_nvm) + 1e-10)
  
  ks_bayesian <- max(abs(actual_cdf - cdf_bayesian))
  ks_em       <- max(abs(actual_cdf - cdf_em))
  ks_nvm      <- max(abs(actual_cdf - cdf_nvm))

  # --- 8. EVALUATION: CHI-SQUARE TEST ---
  # Define number of bins (k).
  n_obs <- length(test_data)
  k <- floor(sqrt(n_obs)) 
  
  # Create bins based on quantiles of test_data to ensure coverage
  breaks <- quantile(test_data, probs = seq(0, 1, length.out = k + 1))
  observed_counts <- as.vector(table(cut(test_data, breaks = breaks, include.lowest = TRUE)))

  # Helper function to get expected counts by integrating the predicted density over bins
  get_expected <- function(pred_y, grid_x, breaks, n_total) {
    # Create a step function/interpolation of the density
    dens_func <- approx(grid_x, pred_y, xout = grid_x, rule = 2)$y
    
    # Calculate probabilities for each bin via trapezoidal integration
    probs <- sapply(1:(length(breaks)-1), function(i) {
      bin_mask <- grid_x >= breaks[i] & grid_x <= breaks[i+1]
      if(sum(bin_mask) < 2) return(1e-10) # Fallback for very narrow bins
      
      # Area under the curve for this bin
      bin_x <- grid_x[bin_mask]
      bin_y <- pred_y[bin_mask]
      width <- diff(bin_x)
      heights <- (bin_y[-1] + bin_y[-length(bin_y)]) / 2
      sum(width * heights)
    })
    
    # Normalize probabilities to sum to 1 and scale to total observations
    probs <- probs / sum(probs)
    return(probs * n_total)
  }

  # --- 9. EVALUATION: BIC (OOS) ---
  # We use the number of parameters (k) from the training phase 
  # and the Log-Likelihood from the test data.
  n_test <- length(test_data)
  
  # k-values
  k_bayesian <- 2  # Mean, Variance
  k_em       <- 2  # G=1: Mean, Variance
  # For NVM, we count components that had significant weight in the training result
  k_nvm      <- if(!is.null(nvm_res)) sum(mix_df[,3] > 1e-4) * 2 else 1

  bic_bayesian_oos <- k_bayesian * log(n_test) - 2 * ll_bayesian
  bic_em_oos       <- k_em * log(n_test) - 2 * ll_em
  bic_nvm_oos      <- k_nvm * log(n_test) - 2 * ll_nvm
  exp_bayesian <- get_expected(bayesian_pred, xj_test, breaks, n_obs)
  exp_em       <- get_expected(em_pred, xj_test, breaks, n_obs)
  exp_nvm      <- get_expected(nvm_pred, xj_test, breaks, n_obs)

  # Calculate Chi-Square Statistic: sum((O - E)^2 / E)
  chisq_bayesian <- sum((observed_counts - exp_bayesian)^2 / (exp_bayesian + 1e-10))
  chisq_em       <- sum((observed_counts - exp_em)^2 / (exp_em + 1e-10))
  chisq_nvm      <- sum((observed_counts - exp_nvm)^2 / (exp_nvm + 1e-10))

  # --- PRINT SUMMARY TABLE ---
  cat(paste("\nOOS ANALYSIS:", name, "\n"))
  cat(sprintf("%-15s | %-12s | %-12s | %-12s\n", "Metric", "Bayesian", "EM", "NVM"))
  cat(rep("-", 60), "\n")
  cat(sprintf("%-15s | %-12.6f | %-12.6f | %-12.6f\n", "MSE (Lower)", mse_bayesian, mse_em, mse_nvm))
  cat(sprintf("%-15s | %-12.2f | %-12.2f | %-12.2f\n", "Log-Lik (Higher)", ll_bayesian, ll_em, ll_nvm))
  cat(sprintf("%-15s | %-12.4f | %-12.4f | %-12.4f\n", "K-S (Lower)", ks_bayesian, ks_em, ks_nvm))
  cat(sprintf("%-15s | %-12.4f | %-12.4f | %-12.4f\n", "Chi-Sq (Lower)", chisq_bayesian, chisq_em, chisq_nvm))
  cat(sprintf("%-15s | %-12.2f | %-12.2f | %-12.2f\n", "BIC (Lower)", bic_bayesian_oos, bic_em_oos, bic_nvm_oos))
  
  mses   <- c(Bayesian = mse_bayesian, EM = mse_em, NVM = mse_nvm)
  chisqs <- c(Bayesian = chisq_bayesian, EM = chisq_em, NVM = chisq_nvm)
  bics   <- c(Bayesian = bic_bayesian_oos, EM = bic_em_oos, NVM = bic_nvm_oos)
  
  cat("OOS Winner (MSE):   ", names(which.min(mses)), "\n")
  cat("OOS Winner (Chi-Sq):", names(which.min(chisqs)), "\n")
  cat("OOS Winner (BIC):   ", names(which.min(bics)), "\n")
  cat("------------------------------------------------------------\n")
}
