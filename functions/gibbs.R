library(MASS)
gibbs_update <- function(
    alpha_sigma, beta_sigma, alpha, 
    w_0, alpha_tau0, beta_tau0, S_inv, 
    tau_i, shape_beta_tau, rate_beta_tau, shape_tau_i,
    w_i, S) {
  
  # Update alpha's variance (sigma2_alpha)
  alpha_sigma_posterior = alpha_sigma + 0.5 * dim(alpha)[1] * dim(alpha)[2]
  beta_sigma_posterior = beta_sigma + 0.5 * sum(alpha^2)
  sigma_alpha_new = rgamma(1, shape = alpha_sigma_posterior, rate = beta_sigma_posterior)
  sigma2_alpha_new = 1.0 / sigma_alpha_new
  
  # Update tau0
  K = length(w_0)
  alpha_tau0_posterior = alpha_tau0 + K / 2
  beta_tau0_posterior = beta_tau0 + 0.5 * (t(w_0) %*% S_inv %*% w_0)
  tau0_new = rgamma(1, shape = alpha_tau0_posterior, rate = beta_tau0_posterior)
  
  # Update beta_tau
  I = length(tau_i)
  shape_beta_tau_posterior = shape_beta_tau + shape_tau_i * I  
  rate_beta_tau_posterior = rate_beta_tau + sum(tau_i)           
  beta_tau_new = rgamma(1, shape = shape_beta_tau_posterior, rate = rate_beta_tau_posterior)
  
  # Update tau_i
  N = ncol(w_i)  # number of regions
  K = nrow(w_i)  # number of spline parameters per region
  tau_i_new = rep(NA, N)
  
  for (i in 1:N) {
    diff = w_i[, i] - w_0
    shape_post = shape_tau_i + K / 2
    rate_post = beta_tau_new+ 0.5 * as.numeric(t(diff) %*% S_inv %*% diff)
    tau_i_new[i] = rgamma(1, shape = shape_post, rate = rate_post)
  }
  
  # Update w0
  total_tau = tau0_new + sum(tau_i_new)
  weighted_matrix = matrix(0, nrow = nrow(w_i), ncol = ncol(w_i)) 
  for (i in 1:length(tau_i_new)) {
    weighted_matrix[, i] = tau_i_new[i] * w_i[, i]
  }
  weighted_sum = rowSums(weighted_matrix)
  mean_posterior = weighted_sum / total_tau
  sigma_posterior = (S / total_tau)
  
  # Sample from posterior distribution of w0
  w0_new = mvrnorm(1, mean_posterior, sigma_posterior)
  
  # Return the updated values
  return(list(
    sigma2_alpha = sigma2_alpha_new,
    tau0 = tau0_new,
    beta_tau= beta_tau_new,
    tau_i = tau_i_new,
    w0 = w0_new
  ))
}


















  