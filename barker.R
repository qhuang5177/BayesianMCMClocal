# A function to evaluate prevalence
# Comment: please note the change in line f_it_vec = 
predict_p = function( beta, w_i, alpha, x, Z,index) {
 
  eta_i            =  Z %*% w_i 

  f_it             = eta_i + alpha  
 
  f_it_vec         = as.vector(f_it)[index]
  
  linear_predictor = f_it_vec + x %*% beta
  
  p_model = 1 / (1 + exp(-linear_predictor))
  
  return(p_model)
}


# function for calculating gradient
# Comment: can find out number of regions from the inputs so removed
# Comment: removed y as input, already in sim_data
# Comment: better to have alpha, beta, and w as inputs alpha is the small fluctuation
# Comment: for alpha, easier to have the variance as input
gradients = function( alpha, beta, w_i, x, Z, sigma2_alpha, S_inv, w0, tau_i, dataset, index,beta_var=10,scale_factor ) {
  
  
  
  
  # Dataset parameters
  nTimes   <- dim(alpha)[1]
  n_region <- dim(alpha)[2]
  K        <- dim(Z)[2]
  
 
    
  
  
  # Current probabilities
  p  = predict_p( beta, w_i, alpha, x, Z, index)
 
  
  
  # Gradient for beta
  # Comment: make 10 a function input, even with a default value
  # Calculating y-p is not expensive but would be ideal to only calculate once
  beta_gradient =scale_factor*(t(x) %*% (dataset$y - p) - (1/beta_var) * beta)
 
  
  # Gradient for alpha
  alpha_it_gradient <- matrix(0, nTimes, n_region)
  for ( i in 1:n_region ) {
    # Comment: calculate only once
    idx1 = dataset$region==i
    for ( t in 1:nTimes) {
      idx                    = (dataset$time==t) & idx1
      grad_sum               = sum( dataset$y[idx]-p[idx] )
      alpha_it_gradient[t,i] = scale_factor*grad_sum - alpha[t,i]/sigma2_alpha
    }
  }
  
  
  # Gradient for w_i
  # Comment: calculate idx1 only once
  # Comment: typo in w0[,i] as there is only one w0:1*6 w_i:6*1
  # Comment: possibly there is a way to avoid triple for loop, this is gonna be slow in R. Just for the future
  w_grad = matrix(0, K, n_region)
 
  for (i in 1:n_region) {
    idx1 = dataset$region==i
    grad_sum = rep(0, K)
    for (t in 1:nTimes) {
      idx = which( (dataset$time==t) & idx1 )
      
      for (j in idx) {
        grad_sum = grad_sum + (dataset$y[j] - p[j]) * Z[t, ]
      }
    }
    w_grad[, i] = scale_factor*grad_sum - tau_i[i] * (S_inv %*% (w_i[, i] - w0))
  }
 
  
  
  # Output
  total_gradient <- c( as.vector(alpha_it_gradient), as.vector(beta_gradient), as.vector(w_grad) )
 
  return(total_gradient)
}



#function for calculating log_posterior  
# Comment: removed y, nTimes, K and n_region from inputs
# Comment: made S_inv an input, no need to invert every time
# Comment: better to directly have sigma2_alpha instead of the SD(varaince not SD)
# Comment: variance of beta better be input
# Comment: w0 is a vector!
# Comment: not sure what the accepted vector was doing
log_posterior = function( alpha, beta, w_i, x, y, Z, sigma2_alpha, S_inv, w0, tau_i, index, beta_var=10) {
  
  # 
  nTimes   = dim(alpha)[1]
  n_region = dim(alpha)[2]
  K        = dim(w_i)[1]
  
 
  # 1. likelihood part
  p_model = predict_p( beta, w_i, alpha, x, Z, index)
  
  #y_model=rbinom(length(p_model),1,p_model)
  loglikelihood = sum(y * log(p_model) + (1 - y) * log(1 - p_model))
  
  # 2. prior part
  
  # prior for alpha
  log_prior_alpha = (- 1/ ( 2*sigma2_alpha) ) *sum(alpha^2)
 
  # prior for beta
  log_prior_beta = - 0.5*sum(beta*beta)/beta_var
  
 
  # prior for w
  log_prior_w_i = 0
  for (j in 1:n_region) {
    diff_j = w_i[,j] - w0
    log_prior_w_i = log_prior_w_i - (tau_i[j]/2) * t(diff_j) %*% S_inv %*% diff_j
  }
  
  # 3. Total log posterior
  logposterior = loglikelihood + log_prior_alpha + log_prior_beta + log_prior_w_i
 
  return(logposterior) 
}   




# Barker update
# Comment: removed current_param
# Comment: n_region can be evaluated from the data
# Comment: removed y from inputs, it's already in sim_data. Fewer inputs, more readability
# Comment: the fitted probabilities are calculated twice, once for likelihood and once for gradient. Could only do one time as it contains expensive exponentials
barker_update = function( beta, alpha, w_i, step_size, x, dataset, Z, sigma2_alpha, S_inv, w0, tau_i, index,batch_index=NULL,scale_factor ) {
  

  # generate noise
  current_param <- c( as.vector(alpha), beta, as.vector(w_i) )
  noise         <- rnorm(length(current_param),step_size,0.1*step_size)
  
 
  
  #select batch data
  x_batch       <- x[batch_index, , drop = FALSE]
  data_batch    <- dataset[batch_index, ]
  index_batch <- index[batch_index,]
  
  
  # current gradient
  beta_x = gradients( alpha, beta, w_i, x_batch, Z, sigma2_alpha, S_inv, w0, tau_i, data_batch, index_batch,beta_var=10,scale_factor)
 
  # Propose new values 
  # Slight speed-up avoiding the loop could be obtained with 2*Bernouli-1 command
  # Comment: typo in the proposal value, b is multiplied by noise, not stepsize
  b    = rep(0,length(current_param))
  prob = rep(0, length(current_param))  
  for (i in 1:length(b)) {
    prob[i] = 1/(1+exp(-beta_x[i]*noise[i]))
    b[i]    = sample(c(-1,1),1,prob = c(1-prob[i],prob[i]))
  }
 
  
  proposal = current_param+b*noise
  
  # Comment: made a functions in utilities to extract the parameters from the vector
  tmp       = from_vector( proposal, dim(alpha)[1], dim(alpha)[2], dim(Z)[2], length(beta) )
  alpha_new = tmp$alpha
  beta_new  = tmp$beta
  w_i_new   = tmp$w_i
  
  
  # 2. Compute log posterior
  log_post_current  = log_posterior( alpha, beta, w_i, x, dataset$y, Z, sigma2_alpha, S_inv, w0, tau_i, index, beta_var=10 ) 
  log_post_proposal = log_posterior( alpha_new, beta_new, w_i_new, x, dataset$y, Z, sigma2_alpha, S_inv, w0, tau_i, index, beta_var=10) 
 
  
  
    
  # 3. Calculate log ratio
  # Comment: could save some calculations by evaluating the formula for the difference
  log_ratio = log_post_proposal - log_post_current
  

  # Gradient at the proposed state
  beta_y=gradients( alpha_new, beta_new, w_i_new, x_batch, Z, sigma2_alpha, S_inv, w0, tau_i, data_batch, index_batch, beta_var=10,scale_factor )
 
  
  # 4. Barker acceptance probability
  # Comment: removed for loop to save some time
  beta1 = -beta_y * (current_param - proposal)
  beta2 = -beta_x * (proposal - current_param)
  correction = sum( -(pmax(beta1,0)+log1p(exp(-abs(beta1)))) + (pmax(beta2,0)+log1p(exp(-abs(beta2)))) )

  
  #  log acceptance ratio
  log_acceptance_ratio = log_ratio + sum(correction)
  # 5. Accept or reject
  u = log(runif(1))
  
  accepted = u < log_acceptance_ratio
 
  if ( accepted ) {
    
  } else {
    beta_new  = beta
    alpha_new = alpha
    w_i_new   = w_i
  }
  
  # Output
  return( list(accepted=accepted, beta=beta_new, alpha=alpha_new, w_i=w_i_new) )
}



