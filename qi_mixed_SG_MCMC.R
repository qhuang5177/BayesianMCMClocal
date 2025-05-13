#Mix SG MCMC

source('qi_Barker method test.R')
source('qi_SG MCMC random.R')

mix_mcmc <- function(nIter, nBurnin, current_param,
                     x, y, Z, sigma_alpha, S, w0, tau_i,
                     nTimes, n_region, K, sim_data,
                     acceptance, batch_idx) {
  
  n_param <- length(current_param)
  posterior_samples <- matrix(NA, nrow = nIter - nBurnin, ncol = n_param)
  
  #initial step size
  step <- rep(0.001, n_param)
  step_avg <- rep(0, n_param)
  
  History <- matrix(NA, nrow = nBurnin, ncol = n_param)
  HistoryAccepted <- matrix(NA, nrow = nBurnin, ncol = n_param)
  
  for (iter in 1:nIter) {
    
    if (iter %% 10 == 0) {
      # using Barker method at the 10th time
      res <- barker_update(current_param, step_size=step, x, y, Z, sigma_alpha, S, w0, tau_i,
                           nTimes, n_region, K, sim_data)
    } else {
      # for the other 9 iterations, use SG
      res <- SG_update(current_param, step, x, y, Z, sigma_alpha, S, w0, tau_i,
                       nTimes, n_region, K, sim_data, batch_idx)
    }
    
    new_param <- res$param
    accepted <- res$accepted
    # record the posterior samples after bunr-in periods
    if (iter > nBurnin) {
      posterior_samples[iter - nBurnin, ] <- new_param
    }
    #adaptive stepsize update under SG MCMC 
    if (iter <= nBurnin && !(iter %% 10 == 0)) {
      History[iter, ] <- step
      HistoryAccepted[iter, ] <- accepted
      
      window <- iter - 1
      if (window >= 1) {
        for (j in 1:n_param) {
          acc_ratio <- mean(HistoryAccepted[1:window, j])
          if (acc_ratio > (acceptance + 0.05)) {
            step[j] <- step[j] * 1.005
          } else if (acc_ratio < (acceptance - 0.05)) {
            step[j] <- step[j] * 0.995
          }
          if (iter > 0.5 * nBurnin) {
            step_avg[j] <- step_avg[j] + 2.0 * step[j] / nBurnin
          }
        }
      }
      
      if (iter == nBurnin) {
        step <- step_avg
      }
    }
    
    current_param <- new_param
  }
  
  return(posterior_samples)
}


#testing
posterior_samples <- mix_mcmc(
  nIter =200,
  nBurnin = 100,
  current_param = init_param,
  x = x,
  y = y,
  Z = Z,
  sigma_alpha = sigma_alpha,
  S = S,
  w0 = w0,
  tau_i = tau_i,
  nTimes = nTimes,
  n_region = n_region,
  K = K,
  sim_data = sim_data,
  acceptance = 0.4,  
  batch_idx = batch_idx
)


