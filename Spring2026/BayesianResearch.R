library(mclust)

################################################################################
# GIBBS SAMPLING FUNCTION
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
# SECTION: TEST INDIVIDUAL RUN
################################################################################
set.seed(123)
test_mu = 100; test_sigma = 5; test_n = 100
y_test = rnorm(test_n, test_mu, test_sigma)
xj_test = seq(min(y_test)-5, max(y_test)+5, length.out = 100)

# 1. Run Models
gibbs_out = run_gibbs(y_test, xj_test)
em_test = densityMclust(y_test, G=1, verbose = FALSE, plot = FALSE)

# 2. Generate Densities
true_dens_test = dnorm(xj_test, test_mu, test_sigma)
bayesian_dens_test = sapply(xj_test, function(x) mean(dnorm(x, gibbs_out$mu, 
                                                            sqrt(gibbs_out$sig2))))
em_dens_test = predict(em_test, xj_test)

# 3. Calculate Performance
mse_bayesian_run = mean((bayesian_dens_test - true_dens_test)^2)
mse_em_run = mean((em_dens_test - true_dens_test)^2)
improvement_run = (mse_em_run - mse_bayesian_run) / mse_em_run * 100

# 4. Print Individual Run Results
cat("--- Individual Run Results ---\n")
cat("MSE Bayesian Model:", mse_bayesian_run, "\n")
cat("MSE EM Algorithm:  ", mse_em_run, "\n")
cat("Bayesian improvement over EM:", round(improvement_run, 2), "%\n\n")

# 5. Diagnostic Plot for Individual Run
par(mfrow=c(1,2)) # Split plot window
# Plot A: Densities
plot(xj_test, true_dens_test, type="l", lwd=2, main="Density Comparison")
lines(xj_test, bayesian_dens_test, col="firebrick", lwd=2)
lines(xj_test, em_dens_test, col="royalblue", lwd=2, lty=2)
legend("topright", legend=c("True", "Bayesian", "EM"), col=c(1, "firebrick", "royalblue"), 
       lty=c(1,1,2), cex=0.7)

# Plot B: Posterior of Sigma
hist(sqrt(gibbs_out$sig2), breaks=30, main="Posterior of Sigma", xlab="Sigma", col="lightgray")
abline(v=test_sigma, col="red", lwd=2) # True Sigma line
par(mfrow=c(1,1)) # Reset plot window


################################################################################
# ITERATIVE COMPARISON (1000 Simulations)

set.seed(123)
n_iterations = 100 # Reduced to 100 for speed, Gibbs is computationally heavier than fixed math
mse_results = data.frame(Bayesian = numeric(n_iterations), EM = numeric(n_iterations))

true_mu = 5; true_sigma = 2; n = 100
xj = seq(true_mu - 4*true_sigma, true_mu + 4*true_sigma, length.out = 100)
true_density = dnorm(xj, mean = true_mu, sd = true_sigma)

for(i in 1:n_iterations) {
  y_sim = rnorm(n, mean = true_mu, sd = true_sigma)
  
  # 1. Gibbs Sampling
  res = run_gibbs(y_sim, n_iter = 1000, burn_in = 200)
  
  # 2. Bayesian Predictive Density
  post_pred = sapply(xj, function(x) mean(dnorm(x, res$mu, sqrt(res$sig2))))
  
  # 3. EM Algorithm
  em_fit = densityMclust(y_sim, G=1, plot = FALSE, verbose = FALSE)
  y_em = predict(em_fit, xj)
  
  # Store MSE
  mse_results$Bayesian[i] = mean((post_pred - true_density)^2)
  mse_results$EM[i] = mean((y_em - true_density)^2)
}

cat("\n--- Final Simulation Results ---\n")
cat("Average MSE Bayesian (Gibbs):", mean(mse_results$Bayesian), "\n")
cat("Average MSE EM:              ", mean(mse_results$EM), "\n")
cat("Bayesian Win Rate:           ", mean(mse_results$Bayesian < mse_results$EM) * 100, "%\n")
