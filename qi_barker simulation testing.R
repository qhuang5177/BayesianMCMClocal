source('qi_simulated data.R')
source('qi_Barker method test.R')

# simulate a dataset for testing
n <- nrow(sim_data)
n_region <- length(unique(sim_data$region))
nTimes <- length(unique(sim_data$time))
K <- 6  


x <- as.matrix(sim_data[, c("x1", "x2", "x3", "x4", "x5", "x6")])
y <- sim_data$y


Z <- matrix(rnorm(nTimes * K), nrow = nTimes, ncol = K)


alpha <- matrix(rnorm(nTimes * n_region), nrow = nTimes, ncol = n_region)
beta <- rnorm(ncol(x))  
w_i <- matrix(rnorm(K * n_region), nrow = K, ncol = n_region)
current_param=c(as.vector(alpha),as.vector(beta),as.vector(w_i))

# Hyperparameter for prior
sigma_alpha <- 1
S <- diag(K)  
w0 <- matrix(0, nrow = K, ncol = n_region)  
tau_i <- rep(1, n_region)  

# Initial value of parameters
params <- c(as.vector(alpha), as.vector(beta), as.vector(w_i))






barker_mcmc <- function(nIter, nBurnin, current_param,
                        x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data,
                        acceptance) {
  n_param <- length(current_param)
  
  
  posterior_samples <- matrix(NA, nrow = nIter - nBurnin, ncol = n_param)
  
  
  step <- rep(0.001, n_param)      # initial step size
 
  # record the history stepsize
  History <- matrix(NA, nrow = nBurnin, ncol = n_param)
  HistoryAccepted <- matrix(NA, nrow = nBurnin, ncol = n_param)
  
  # record average step size
  step_avg <- rep(0, n_param)
  
  for (iter in 1:nIter) {
    
    res <- barker_update(current_param, step, x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data)
    new_param <- res$param
    accepted <- res$accepted
    
    # After Burn-in save the posterior parameter
    if (iter > nBurnin) {
      posterior_samples[iter - nBurnin, ] <- new_param
    }
    
    # adaptive step size 
    if (iter <= nBurnin) {
      History[iter, ] <- step
      HistoryAccepted[iter, ] <- accepted * 1  # （TRUE/FALSE → 1/0）
      
      # stepsize revise
      window <- iter - 1
      if (window >= 1) {
        for (j in 1:n_param) {
          acceptance_ratio <- mean(HistoryAccepted[1:window, j])
          
          if (acceptance_ratio > acceptance + 0.05) {
            step[j] <- step[j] * 1.005
          } else if (acceptance_ratio < acceptance - 0.05) {
            step[j] <- step[j] * 0.995
          }
          
          # update average
          if (iter > 0.5 * nBurnin) {
            step_avg[j] <- step_avg[j] + 2.0 * step[j] / nBurnin
          }
        }
      }
      
      # last step of burn-in，replace it with average step size
      if (iter == nBurnin) {
        step <- step_avg
      }
    }
    
   
    current_param <- new_param
  }
  
  return(posterior_samples)
}











#Testing code



samples=barker_mcmc(2000, 1000, params,
                    x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data,
                    acceptance = 0.4)






alpha_it_len = nTimes * n_region
beta_len = ncol(x)  
w_i_len = n_region * K



# extract all samples for alpha_it 
alpha_chain = samples[, 1:alpha_it_len]

#extract all samples for beta 
beta_chain = samples[, (alpha_it_len + 1):(alpha_it_len + beta_len)]

# extract all samples for  w_i 
w_i_chain = samples[, (alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]



#trace plot 
plot(alpha_chain[,1], type = "l")
plot(beta_chain[,1], type = "l")
plot(w_i_chain[,1], type = "l")





