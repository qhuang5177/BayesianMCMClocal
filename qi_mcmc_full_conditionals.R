library(MASS) 
update_sigma_random_effects = function( alpha_sigma, beta_sigma, alpha )#alpha corresponding to alpha_it
{
  # Calculate the parameters of the full conditional
  alpha_sigma_posterior = alpha_sigma + 0.5*dim(alpha)[1]*dim(alpha)[2]
  beta_sigma_posterior = beta_sigma + 0.5*sum(alpha^2)
  # Draw the new value
  sigma_alpha_new = rgamma( 1, shape=alpha_sigma_posterior, rate=beta_sigma_posterior )
  return( 1.0/sigma_alpha_new )
  
}


update_tau0 = function(w_0, alpha_tau0,beta_tau0) {
  #calculate the parameters of the full conditional
  K = length(w_0)
  alpha_tau0_posterior = alpha_tau0+ K / 2
  beta_tau0_posterior = beta_tau0+ 0.5 * as.numeric(t(w_0) %*% solve(S) %*% w_0)
  # Draw the new value
  tau0_new = rgamma(1, shape = alpha_tau0_posterior, rate = beta_tau0_posterior)
  return(tau0_new)
}

update_beta_tau <- function(tau_i, shape_beta_tau, rate_beta_tau , shape_tau_i ) {
  I=length(tau_i)
  # calculate the parameters of the full conditional
  shape_beta_tau_posterior = shape_beta_tau + shape_tau_i * I  
  rate_beta_tau_posterior = rate_beta_tau  + sum(tau_i)           
  # Draw the new value
  beta_tau_new = rgamma(1, shape = shape_beta_tau_posterior, rate = rate_beta_tau_posterior)
  return(beta_tau_new)
}

update_tau_i <- function(w_i, w_0, shape_tau_i, tau_i_rate) {
    K =length(w_i)
   # calculate the parameters of the full conditional
   tau_i_shape_posterior<-shape_tau_i+K/2
   tau_i_rate_posterior <-tau_i_rate+0.5* as.numeric(t(w_i - w_0) %*% solve(S) %*% (w_i - w_0))
   # Draw the new value
   tau_i_new=rgamma(1, shape = tau_i_shape_posterior, rate = tau_i_rate_posterior)
   return(tau_i_new)
}


update_w_0 <- function(tau_0, tau_i, w_i, S) {
  total_tau <- tau_0 + sum(tau_i)
  weighted_list <- mapply(function(w, tau) tau * w, w_i, tau_i, SIMPLIFY = FALSE)
  weighted_sum <- Reduce("+", weighted_list)
   #Posterior mean 
  mean_posterior <- weighted_sum / total_tau
  #  Posterior sigma 
  sigma_posterior <- (S / total_tau)
  #  Sample from posterior
  w_0_new <- MASS::mvrnorm(1, mean_posterior, sigma_posterior)
  
  return(w_0_new)
}







