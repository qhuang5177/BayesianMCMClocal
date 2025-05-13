nTimes=14
K=6
n_region=23





#function for calculating prevelance
predict_p = function(beta, w_i,alpha,x, Z) {
  
  eta_i =  Z %*% w_i  
  f_it = eta_i + alpha          
  f_it_vec = as.vector(f_it)
  linear_predictor = f_it_vec + x %*% beta
  
  p_model = 1 / (1 + exp(-linear_predictor))
  
  return(p_model)
}


#function for calculating gradient
gradients = function(params,x, y,Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K,sim_data) {
  
  alpha_it_len = nTimes * n_region
  beta_len = ncol(x)  
  w_i_len = n_region * K
  
  #The input may be a vector, now reshape it to a matrix for calculating gradient
  alpha = params[1:alpha_it_len]
  alpha = matrix(alpha, nrow = nTimes, ncol = n_region)
  
  beta = params[(alpha_it_len + 1):(alpha_it_len + beta_len)]
  
  w_i = params[(alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]
  w_i = matrix(w_i, nrow = K, ncol = n_region)
  
  p  = predict_p(beta,w_i,alpha,x, Z)
  S_inv = solve(S)
  
  #Gradient for alpha_it
  alpha_it_gradient <- matrix(0, nTimes, n_region)
  
  for (i in 1:n_region) {
    for (t in 1:nTimes) {
      idx = which(sim_data$region == i & sim_data$time == t)
      grad_sum = sum(y[idx] - p[idx])
      alpha_it_gradient[t, i] = grad_sum - alpha[t, i] / (sigma_alpha^2)
    }
  }
  
  # Gradient for beta
  beta_gradient = t(x) %*% (y - p) - (1/10) * beta
  
  # Gradient for w_i
  w_grad = matrix(0, K, n_region)

  for (i in 1:n_region) {
    grad_sum = rep(0, K)
    for (t in 1:nTimes) {
      idx = which(sim_data$region == i & sim_data$time == t)  
      for (j in idx) {
        grad_sum = grad_sum + (y[j] - p[j]) * Z[t, ]
      }
    }
    w_grad[, i] = grad_sum - tau_i[i] * (S_inv %*% (w_i[, i] - w0[,i]))
  }
  
  total_gradient <- c(as.vector(alpha_it_gradient), as.vector(beta_gradient), as.vector(w_grad))
  
  return(total_gradient)
}
  
#function for calculating log_posterior  
log_posterior = function(params,x,y, Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K) {
  
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
  
  p_model = predict_p(beta,w_i,alpha,x, Z)
  #y_model=rbinom(length(p_model),1,p_model)
  
  loglikelihood = sum(y * log(p_model) + (1 - y) * log(1 - p_model))
  
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
  logposterior = loglikelihood + log_prior_alpha + log_prior_beta + log_prior_w_i
  
  return(as.numeric(logposterior)) 
  
}   
  


update_stepsize <- function(mcmc_idx, nBurn, target_lower, target_upper,
                            accept_matrix, stepsize) {
  window <- nrow(accept_matrix) - 1  # length of windows
  nVars <- ncol(accept_matrix)      # length of variable
  
  for (j in 1:nVars) {
    # 1. calcualte the acceptance ratio
    acceptance_ratio <- mean(accept_matrix[1:window, j])
    
    # 2. adjust stepszie
    if (acceptance_ratio > target_upper) {
      stepsize[j] <- stepsize[j] * 1.005
    } else if (acceptance_ratio < target_lower) {
      stepsize[j] <- stepsize[j] * 0.995
    }
    
    # 3. calculate the average（only for the half of brun-in period）
    if (mcmc_idx > 0.5 * nBurn) {
      average[j] <- average[j] + 2.0 * stepsize[j] / nBurn
    }
  }
  
  # 4. final iteration in burn-in 
  if (mcmc_idx == nBurn) {
    stepsize <- average
  }
  
  return(list(stepsize = stepsize, average = average))
}



  
  
  
  
#calculation  barker method
barker_update = function(current_param,step_size,x, y, Z, sigma_alpha, S, w0, tau_i, nTimes, n_region, K, sim_data) {
  
  
 #generate noise
  noise<-rnorm(length(current_param),0,step_size)
  #gradient
  b=rep(0,length(current_param))
  prob = rep(0, length(current_param))  
  beta_x=gradients(current_param,x,y,Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K,sim_data)
 
  for (i in 1:length(b)) {
    prob[i]=1/(1+exp(-beta_x[i]*noise[i]))
    b[i]=sample(c(-1,1),1,prob = c(1-prob[i],prob[i]))
  }
  
  
  
  #proposal
  proposal=current_param+b*step_size
  
  # 2. Compute log posterior
  log_post_current = log_posterior(current_param,x, y,Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K)
  log_post_proposal = log_posterior(proposal,x,y, Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K)
  
  # 3. Calculate log ratio
  log_ratio = log_post_proposal - log_post_current
  
  # 4. Barker acceptance probability
  correction_logterms = rep(0, length(current_param))
  
  beta_y=gradients(proposal,x, y,Z, sigma_alpha, S, w0,tau_i,nTimes, n_region, K,sim_data)
  

    for (i in 1:length(current_param)) {
      beta1 = -beta_y[i] * (current_param[i] - proposal[i])
      beta2 = -beta_x[i] * (proposal[i] - current_param[i])
      
      # Barker trick: log(1 / (1 + exp(-β))) = -log(1 + exp(-β))
      log_term1 = -(pmax(beta1, 0) + log1p(exp(-abs(beta1))))
      log_term2 = +(pmax(beta2, 0) + log1p(exp(-abs(beta2))))
      
      correction_logterms[i] = log_term1 + log_term2
    }
  
  #  log acceptance ratio
  log_acceptance_ratio = log_ratio + sum(correction_logterms)
  # 5. Accept or reject
  u = log(runif(1))
  accepted = u < log_acceptance_ratio
 
  
  accepted_vector <- as.numeric(b != 0 & accepted)
  
  next_point <- if (accepted) proposal else current_param
  
  
  return(list(param = next_point, accepted = accepted_vector))
}



