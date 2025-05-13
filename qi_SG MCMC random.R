nTimes=14
K=6
n_region=23






#function for calculating prevelance
SG_predict_p = function(beta, w_i, alpha, x_batch, Z, sim_data) {
  region_batch <- sim_data$region[batch_idx]
  time_batch <- sim_data$time[batch_idx]
  
  eta_mat = Z %*% w_i       
  f_mat = eta_mat + alpha   # (nTimes x nRegion)
  
  idx = cbind(time_batch, region_batch)
  f_batch = f_mat[idx]
  
  linear_predictor = f_batch + x_batch %*% beta
  p_model = 1 / (1 + exp(-linear_predictor))
  
  return(p_model)
}

#function for calculating gradient
SG_gradients = function(params, x, y, Z, sigma_alpha, S, w0, tau_i,
                     nTimes, n_region, K, sim_data, batch_idx) {
  
  # batch 
  x_batch = x[batch_idx, ]
  y_batch = y[batch_idx]
  sim_data_batch = sim_data[batch_idx, ]
  region_batch <- sim_data$region[batch_idx]
  time_batch <- sim_data$time[batch_idx]
  
  
  alpha_it_len = nTimes * n_region
  beta_len = ncol(x)
  w_i_len = n_region * K
  
  alpha = matrix(params[1:alpha_it_len], nrow = nTimes, ncol = n_region)
  beta = params[(alpha_it_len + 1):(alpha_it_len + beta_len)]
  w_i = matrix(params[(alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)], nrow = K, ncol = n_region)
  
  p = SG_predict_p(beta, w_i, alpha, x_batch, Z,sim_data)
  S_inv = solve(S)
  
  # Gradient for alpha
  alpha_it_gradient <- matrix(0, nTimes, n_region)
  for (i in 1:n_region) {
    for (t in 1:nTimes) {
      idx = which(sim_data_batch$region == i & sim_data_batch$time == t)
      grad_sum = sum(y_batch[idx] - p[idx])
      alpha_it_gradient[t, i] = grad_sum - alpha[t, i] / (sigma_alpha^2)
    }
  }
  
  # Gradient for beta
  beta_gradient = t(x_batch) %*% (y_batch - p) - (1/10) * beta
  
  # Gradient for w_i
  w_grad = matrix(0, K, n_region)
  for (i in 1:n_region) {
    grad_sum = rep(0, K)
    for (t in 1:nTimes) {
      idx = which(sim_data_batch$region == i & sim_data_batch$time == t)
      for (j in idx) {
        grad_sum = grad_sum + (y_batch[j] - p[j]) * Z[t, ]
      }
    }
    w_grad[, i] = grad_sum - tau_i[i] * (S_inv %*% (w_i[, i] - w0[, i]))
  }
  scaling_factor = nrow(x) / length(batch_idx)  # 即 N / B
  
  # 然后在 return 前统一加：
  total_gradient <- scaling_factor * c(as.vector(alpha_it_gradient), 
                                       as.vector(beta_gradient), 
                                       as.vector(w_grad))
  
  return(total_gradient)
}

  
 
  
  
  
  
  
  
  
  
  

#function for calculating log_posterior  
SG_log_posterior = function(params,x,y, Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K,batch_idx) {
  
  
  x_batch = x[batch_idx, ]
  y_batch = y[batch_idx]
 
  
  ##The input may be a vector, now reshape it to a matrix 
  alpha_it_len = nTimes * n_region
  beta_len = ncol(x)  
  w_i_len = n_region * K
  
  
  alpha = params[1:alpha_it_len]
  alpha = matrix(alpha, nrow = nTimes, ncol = n_region)
  
  beta = params[(alpha_it_len + 1):(alpha_it_len + beta_len)]
  
  w_i = params[(alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]
  w_i = matrix(w_i, nrow = K, ncol = n_region)
  
  S_inv = solve(S)
  
  
  
  # 1. likelihood part
  
  p_model = SG_predict_p(beta,w_i,alpha,x_batch, Z,sim_data)
  #y_model=rbinom(length(p_model),1,p_model)
  
  loglikelihood = sum(y_batch * log(p_model) + (1 - y_batch) * log(1 - p_model))
  scaling_factor = nrow(x) / length(batch_idx) 
  total_loglikelihood=scaling_factor*loglikelihood
  # 2. prior part
  
  # prior for alpha
  log_prior_alpha = (- 1/ ( 2*sigma_alpha^2)) *sum(alpha^2)
  
  
  # prior for beta
  log_prior_beta = - sum(beta*beta)/20
  
  # prior for w
  log_prior_w_i = 0
  for (j in 1:n_region) {
    diff_j = w_i[,j] - w0[,j]
    log_prior_w_i = log_prior_w_i - (tau_i[j]/2) * t(diff_j) %*% S_inv %*% diff_j
  }
  
  
  
  
  # 3. Total log posterior
  logposterior = total_loglikelihood + log_prior_alpha + log_prior_beta + log_prior_w_i
  
  return(as.numeric(logposterior)) 
  
}   










#calculation  barker method
SG_update = function(current_param, step_size, x, y, Z, sigma_alpha, S, w0, tau_i,
                         nTimes, n_region, K, sim_data, batch_idx) {
  
  noise <- rnorm(length(current_param), 0, step_size)
  b = rep(0, length(current_param))
  prob = rep(0, length(current_param))
  
  beta_x = SG_gradients(current_param, x, y, Z, sigma_alpha, S, w0, tau_i,
                     nTimes, n_region, K, sim_data, batch_idx)
  
  for (i in 1:length(b)) {
    prob[i] = 1 / (1 + exp(-beta_x[i] * noise[i]))
    b[i] = sample(c(-1, 1), 1, prob = c(1 - prob[i], prob[i]))
  }
  
  proposal = current_param + b * step_size
  
  log_post_current = SG_log_posterior(current_param, x, y, Z, sigma_alpha, S, w0, tau_i,
                                   nTimes, n_region, K,batch_idx) 
  log_post_proposal = SG_log_posterior(proposal, x, y, Z, sigma_alpha, S, w0, tau_i,
                                    nTimes, n_region, K,batch_idx)
  
  log_ratio = log_post_proposal - log_post_current
  
  correction_logterms = rep(0, length(current_param))
  beta_y = SG_gradients(proposal, x, y, Z, sigma_alpha, S, w0, tau_i,
                     nTimes, n_region, K, sim_data, batch_idx)
  
  for (i in 1:length(current_param)) {
    beta1 = -beta_y[i] * (current_param[i] - proposal[i])
    beta2 = -beta_x[i] * (proposal[i] - current_param[i])
    log_term1 = -(pmax(beta1, 0) + log1p(exp(-abs(beta1))))
    log_term2 = +(pmax(beta2, 0) + log1p(exp(-abs(beta2))))
    correction_logterms[i] = log_term1 + log_term2
  }
  
  log_acceptance_ratio = log_ratio + sum(correction_logterms)
  u = log(runif(1))
  accepted <- u < log_acceptance_ratio
  
  next_point <- if (accepted) proposal else current_param
  accepted_vector <- as.numeric(b != 0 & accepted)
  
  return(list(param = next_point, accepted = accepted_vector))
}






