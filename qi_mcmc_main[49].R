stoch_gradient_logistic = function(nIters=12000, nThin=10, initial_values, alpha_sigma, beta_sigma,alpha_tau0,beta_tau0,shape_beta_tau,rate_beta_tau,shape_tau_i,tau_i_rate )
{
  
  
  # 
  source('qi_mcmc_full_conditionals.R')
  #nT = length( table(Y$time) )
  #N = length( table(Y$region) )
  
  # Storage
  sigma_alpha_mcmc = rep(NA,nIters/nThin)
  tau_0_mcmc=rep(NA,nIters/nThin)
  beta_tau_mcmc=rep(NA,nIters/nThin)
  tau_i_mcmc = matrix(NA,nrow = nIters / nThin, ncol = length(initial_values$tau_i))
  w_0_mcmc = matrix(NA, nrow = nIters / nThin, ncol = length(initial_values$w_0))
  
  # Set initial values
  alpha = initial_values$alpha
  w_0 = initial_values$w_0
  w_i = initial_values$w_i
  tau_i = initial_values$tau_i
  tau_0 = initial_values$tau_0

  
  nT = 14
  N = 23

  idx = 1
  for ( i in 1:nIters ) {
     
    # Update variance of random effects 
    #sigma_alpha = update_sigma_random_effects( alpha_sigma, beta_sigma, alpha )
    
    #Update precision of spline term tao_0
    #tau_0=update_tau0(w_0, alpha_tau0,beta_tau0,S)
    
    #update beta_tau
    #beta_tau=update_beta_tau (tau_i, shape_beta_tau, rate_beta_tau , shape_tau_i ) 
    
    #update tau_i
    #tau_i=update_tau_i(w_i, w_0, shape_tau_i, tau_i_rate)
    # Update omega_zero
    w_0=update_w_0(tau_0,tau_i,w_i,S)
  
    
    # save when appropriate
    if ( i%%nThin == 0 ){
      sigma_alpha_mcmc[idx] = sigma_alpha#output from function update_sigma_random _effect
      tau_0_mcmc[idx]=tau_0
      beta_tau_mcmc[idx]=beta_tau
      tau_i_mcmc[idx,]=tau_i
      w_0_mcmc[idx,]=w_0
      idx = idx + 1
    }
  }
  
  # Output 
  return(
    list( sigma_alpha_mcmc=sigma_alpha_mcmc,tau_0_mcmc=tau_0_mcmc,beta_tau_mcmc=beta_tau_mcmc,tau_i_mcmc=tau_i_mcmc,w_0_mcmc=w_0_mcmc)
  )
}

