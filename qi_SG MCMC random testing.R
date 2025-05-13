source('qi_SG MCMC random.R')
source('qi_simulated data.R')

#simulation dataset for testing 
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












SG_mcmc <- function(nIter, nBurnin, current_param,
                        x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data,
                        acceptance,batch_idx) {
  n_param <- length(current_param)
  
  
  posterior_samples <- matrix(NA, nrow = nIter - nBurnin, ncol = n_param)
  
  # initial parameters
  step <- rep(0.001, n_param)      
  
  
  # record the pace for history 
  History <- matrix(NA, nrow = nBurnin, ncol = n_param)
  HistoryAccepted <- matrix(NA, nrow = nBurnin, ncol = n_param)
  
  # record the average step size 
  step_avg <- rep(0, n_param)
  
  for (iter in 1:nIter) {
    
    res <- SG_update(current_param, step, x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data,batch_idx)
    new_param <- res$param
    accepted <- res$accepted
    
    # After Burn-in 
    if (iter > nBurnin) {
      posterior_samples[iter - nBurnin, ] <- new_param
    }
    
    # adaptive step size
    if (iter <= nBurnin) {
      History[iter, ] <- step
      HistoryAccepted[iter, ] <- accepted * 1  # （TRUE/FALSE → 1/0）
      
      # 
      window <- iter - 1
      if (window >= 1) {
        for (j in 1:n_param) {
          acceptance_ratio <- mean(HistoryAccepted[1:window, j])
          
          if (acceptance_ratio > acceptance + 0.05) {
            step[j] <- step[j] * 1.005
          } else if (acceptance_ratio < acceptance - 0.05) {
            step[j] <- step[j] * 0.995
          }
          
          # update step_avg
          if (iter > 0.5 * nBurnin) {
            step_avg[j] <- step_avg[j] + 2.0 * step[j] / nBurnin
          }
        }
      }
      
      # last step in burn-in 
      if (iter == nBurnin) {
        step <- step_avg
      }
    }
    
  
    current_param <- new_param
  }
  
  return(posterior_samples)
}


#testing the code 
set.seed(123)

# simulate batch index
n <- nrow(sim_data)
batch_size <- 100
batch_idx <- sample(1:n, batch_size)


# initial parameters
params <- c(as.vector(alpha), as.vector(beta), as.vector(w_i))









posterior_samples <- SG_mcmc(
  nIter =2000,
  nBurnin = 1000,
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





alpha_it_len = nTimes * n_region
beta_len = ncol(x)  
w_i_len = n_region * K



# extract all the posterior samples for  alpha_it 
alpha_chain = posterior_samples[, 1:alpha_it_len]

# extract all the posterior samples for  beta 
beta_chain = posterior_samples[, (alpha_it_len + 1):(alpha_it_len + beta_len)]

# extract all the posterior samples for w_i
w_i_chain = posterior_samples[, (alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]


#plot the trace plot 
plot(alpha_chain[,1], type = "l")

plot(beta_chain[,1], type = "l")
plot(w_i_chain[,1], type = "l")









