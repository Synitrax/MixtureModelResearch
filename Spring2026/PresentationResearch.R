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
    phi_full = phi_matrix, # The matrix you requested
    mu = mu_final,
    sig2 = sig2_final
  ))
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
dev.new() # Opens one clean graphics window
par(mfrow = c(2, 3))

# Use the full range for both training and testing
for(name in names(stock_list)) {
  # Full data for both
  train_data <- stock_list[[name]]
  test_data  <- train_data # Testing on what we trained on
  
  # 1. Define Evaluation Grid based on the full data range
  xj_test = seq(min(test_data), max(test_data), length.out = 1000)
  
  # --- MODEL 1: Bayesian (Gibbs) ---
  # Averaging over the posterior samples
  gibbs_train = run_gibbs(train_data, xj_test)
  bayesian_pred = gibbs_train$y
  
  # --- MODEL 2: EM (Mclust) ---
  # Single Gaussian fit via EM
  em_train = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  em_pred  = predict(em_train, xj_test)
  
  # --- MODEL 3: NVMunmix ---
  mu_grid <- seq(min(train_data), max(train_data), length.out = 25)
  r_grid  <- seq(sd(train_data)*0.3, sd(train_data)*1.2, length.out = 15)
  
  nvm_res <- tryCatch({
    NVMunmix(train_data, b=0.4, rgrid=r_grid, mugrid=mu_grid, tol=0.001)
  }, error = function(e) return(NULL))
  
  if(!is.null(nvm_res)) {
    mix_df <- nvm_res[[1]]
    nvm_pred <- sapply(xj_test, function(x) {
      sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2]))
    })
  } else {
    nvm_pred <- rep(0, length(xj_test))
  }
  # 5. Plotting (Add this inside the loop)
  plot(density(test_data), main=paste("In-Sample Fit:", name), lwd=1, col="gray")
  lines(xj_test, bayesian_pred, col="firebrick", lwd=2) # Bayesian
  lines(xj_test, em_pred, col="royalblue", lwd=2, lty=2) # EM
  lines(xj_test, nvm_pred, col="darkgreen", lwd=2, lty=3) # NVM

  # Optional: Refresh legend on one of the plots
  legend("topright", legend=c("Empirical", "Bayesian", "EM", "NVM"), 
         col=c("gray", "firebrick", "royalblue", "darkgreen"), lty=c(1,1,2,3), cex=0.7)
                      
  # --- 3. EVALUATE ---
  # Generate true density of the same data for comparison
  # IMPORTANT: Set n=200 to match xj_test length for a fair MSE
  true_dens_in = density(test_data, from=min(xj_test), to=max(xj_test), n=1000)$y
  
  mse_bayesian = mean((bayesian_pred - true_dens_in)^2)
  mse_em       = mean((em_pred - true_dens_in)^2)
  mse_nvm      = mean((nvm_pred - true_dens_in)^2)
  
  mses <- c(Bayesian = mse_bayesian, EM = mse_em, NVM = mse_nvm)
  winner <- names(which.min(mses))
  
  cat(paste("Stock:", name, "\n"))
  cat("In-Sample MSE Bayesian:", round(mse_bayesian, 6), "\n")
  cat("In-Sample MSE EM:      ", round(mse_em, 6), "\n")
  cat("In-Sample MSE NVM:     ", round(mse_nvm, 6), "\n")
  cat("In-Sample Winner:      ", winner, "\n")
  cat("----------------------------\n")

  # 1. Create Bins (Frequency) from the test data
  # 'plot=FALSE' ensures we don't open more windows
  h <- hist(test_data, breaks = "FD", plot = FALSE)
  observed <- h$counts
  bin_breaks <- h$breaks

  # 2. Function to get Expected Counts for a model
  get_expected <- function(pred_func_or_values, breaks, total_n) {
    # We integrate the density over each bin interval
    # For simplicity, we can use the midpoint approximation
    midpoints <- breaks[-1] - diff(breaks)/2
  
    # If using the pred vectors we already made:
    # We need to map the pred values to the specific bin midpoints
    # This is a simplified version for your presentation:
    bin_probs <- sapply(1:(length(breaks)-1), function(i) {
      # Use the area of the bin: height * width
      # We'll use the mean density in that bin range
      mean(pred_func_or_values[xj_test >= breaks[i] & xj_test <= breaks[i+1]]) * diff(breaks)[i]
    })
  
    # Ensure probabilities don't result in 0 (to avoid division by zero)
    bin_probs[bin_probs <= 0] <- 1e-10
    return(bin_probs * total_n)
  }
  # 3. Calculate Chi-Square for each
  exp_bayesian <- get_expected(bayesian_pred, bin_breaks, length(test_data))
  exp_em       <- get_expected(em_pred, bin_breaks, length(test_data))
  exp_nvm      <- get_expected(nvm_pred, bin_breaks, length(test_data))

  chisq_bayesian <- sum((observed - exp_bayesian)^2 / exp_bayesian)
  chisq_em       <- sum((observed - exp_em)^2 / exp_em)
  chisq_nvm      <- sum((observed - exp_nvm)^2 / exp_nvm)

  # --- Print Results ---
  cat("Chi-Square Results:\n")
  cat("Bayesian:", round(chisq_bayesian, 2), "\n")
  cat("EM:      ", round(chisq_em, 2), "\n")
  cat("NVM:     ", round(chisq_nvm, 2), "\n")
}


############################################################################
dev.new() # Opens one clean graphics window
par(mfrow = c(2, 3))

split_date <- as.Date("2025-04-03")

cat("--- Out-of-Sample Performance: Bayesian vs EM vs NVM ---\n\n")

for(name in names(stock_list)) {
  # Extract numeric vector and match with dates
  full_series <- stock_list[[name]]
  
  # Split based on the dates vector
  train_data <- full_series[returns_dates < split_date]
  test_data  <- full_series[returns_dates >= split_date]
  
  # Define Evaluation Grid based on TEST data range
  xj_test = seq(min(test_data), max(test_data), length.out = 1000)
  
  # --- MODEL 1: Bayesian (Gibbs) ---
  gibbs_train = run_gibbs(train_data, xj_test)
  bayesian_pred = gibbs_train$y
  
  # --- MODEL 2: EM (Mclust) ---
  em_train = densityMclust(train_data, G=1, verbose = FALSE, plot = FALSE)
  em_pred  = predict(em_train, xj_test)
  
  # --- MODEL 3: NVMunmix ---
  # We define the mu and sigma grids dynamically based on the stock's profile
  mu_grid <- seq(min(train_data), max(train_data), length.out = 20)
  r_grid  <- seq(sd(train_data)*0.5, sd(train_data)*1.5, length.out = 10)
  
  nvm_res <- tryCatch({
    # Call your original NVMunmix function
    NVMunmix(train_data, b=0.4, rgrid=r_grid, mugrid=mu_grid, tol=0.001)
  }, error = function(e) return(NULL))
  
  # Calculate NVM Density if successful
  if(!is.null(nvm_res)) {
    mix_df <- nvm_res[[1]]
    nvm_pred <- sapply(xj_test, function(x) {
      sum(mix_df[,3] * dnorm(x, mix_df[,1], mix_df[,2]))
    })
  } else {
    nvm_pred <- rep(0, length(xj_test))
  }

  # 5. Plotting (Add this inside the loop)
  plot(density(test_data), main=paste("OOS Predict:", name), lwd=1, col="gray")
  lines(xj_test, bayesian_pred, col="firebrick", lwd=2) # Bayesian
  lines(xj_test, em_pred, col="royalblue", lwd=2, lty=2) # EM
  lines(xj_test, nvm_pred, col="darkgreen", lwd=2, lty=3) # NVM

  # Optional: Refresh legend on one of the plots
  legend("topright", legend=c("Empirical", "Bayesian", "EM", "NVM"), 
         col=c("gray", "firebrick", "royalblue", "darkgreen"), lty=c(1,1,2,3), cex=0.7)

  # --- 3. EVALUATE ---
  true_dens_future = density(test_data, from=min(xj_test), to=max(xj_test), n=1000)$y
  
  mse_bayesian = mean((bayesian_pred - true_dens_future)^2)
  mse_em       = mean((em_pred - true_dens_future)^2)
  mse_nvm      = mean((nvm_pred - true_dens_future)^2)
  
  mses <- c(Bayesian = mse_bayesian, EM = mse_em, NVM = mse_nvm)
  winner <- names(which.min(mses))
  
  cat(paste("Stock:", name, "\n"))
  cat("MSE Bayesian:", round(mse_bayesian, 6), "\n")
  cat("MSE EM:      ", round(mse_em, 6), "\n")
  cat("MSE NVM:     ", round(mse_nvm, 6), "\n")
  cat("Winner:      ", winner, "\n")
  cat("----------------------------\n")

  h <- hist(test_data, breaks = "FD", plot = FALSE)
  observed <- h$counts
  bin_breaks <- h$breaks

  # 2. Function to get Expected Counts for a model
  get_expected <- function(pred_func_or_values, breaks, total_n) {
    # We integrate the density over each bin interval
    # For simplicity, we can use the midpoint approximation
    midpoints <- breaks[-1] - diff(breaks)/2
  
    # If using the pred vectors we already made:
    # We need to map the pred values to the specific bin midpoints
    # This is a simplified version for your presentation:
    bin_probs <- sapply(1:(length(breaks)-1), function(i) {
      # Use the area of the bin: height * width
      # We'll use the mean density in that bin range
      mean(pred_func_or_values[xj_test >= breaks[i] & xj_test <= breaks[i+1]]) * diff(breaks)[i]
    })
  
    # Ensure probabilities don't result in 0 (to avoid division by zero)
    bin_probs[bin_probs <= 0] <- 1e-10
    return(bin_probs * total_n)
  }
  # 3. Calculate Chi-Square for each
  exp_bayesian <- get_expected(bayesian_pred, bin_breaks, length(test_data))
  exp_em       <- get_expected(em_pred, bin_breaks, length(test_data))
  exp_nvm      <- get_expected(nvm_pred, bin_breaks, length(test_data))

  chisq_bayesian <- sum((observed - exp_bayesian)^2 / exp_bayesian)
  chisq_em       <- sum((observed - exp_em)^2 / exp_em)
  chisq_nvm      <- sum((observed - exp_nvm)^2 / exp_nvm)

  # --- Print Results ---
  cat("Chi-Square Results:\n")
  cat("Bayesian:", round(chisq_bayesian, 2), "\n")
  cat("EM:      ", round(chisq_em, 2), "\n")
  cat("NVM:     ", round(chisq_nvm, 2), "\n")
}
