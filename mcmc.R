
# 
# Comment: for as many inputs as you can, good to set default values in case the user forgets
# Comment: I have removed current_parameter as an input, I recommend making the function more general, for all parameters not just ones updated with Barker
# Comment: Number of regions does not need to be an input, it can be calculated
# Comment: x and y does not need to be an input, right? Already in sim_data
# Comment: currently code only works if regions are coded 1-N and times 1-T. Could make it more general in the future
# Comment: for adaptation window, better to look at 50-200 previous iterations
barker_mcmc <- function(
    nIter, 
    nBurnin, 
    init=NULL,  # Comment: good to give option to provide initial values. Default is NULL, need to make a function that generates some values
    nKnots = 6,
    dataset,                  # This is the dataset the contains region, time and covariates
    target_acceptance=0.55,   # Comment: slightly changed the name of this variable
    covariates,               # Which columns have the covariates 
    nThin=1,
    alpha_sigma = 0.01,         
    beta_sigma = 0.01,
    alpha_tau0 = 10, 
    beta_tau0 =  1,
    shape_beta_tau = 1, 
    rate_beta_tau = 100,
    shape_tau_i = 10,         # Thinning: if there are too many MCMC iterations, you can't save all of the samples
    window=100, 
    batch_type = "full",      # Type of gradient: full, stratified, stratified_min
    hybrid = FALSE,           # Add some regular Barker updates every now and then   
    barker_frequency = 10,    # How often to do Barker in hybrid
    prop = 0.2,               # Proportion of sample to be drawn if using batch
    min_per_group = 10        # Minimum samples per group for stratified sampling
    )
{
  
  # Some global variables to be used 
  n_region <- length( table(dataset$region) )
  nTimes   <- length( table(dataset$time) )
  x        <- as.matrix(dataset[,c(covariates)])
  y        <- dataset$y
  region   <- dataset$region
  time     <- dataset$time
  P        <- dim(x)[2]
  
  
  # To help with the evaluation of gradients
  index <- individuals_index( dataset$region, dataset$time )
 
    
  # Set up the spline here
  tmp <- spline_setup( nTimes, nKnots )
  Z   <- tmp$Z
  # Comment: No need to invert many time within the function. You can invert here once
  S_inv <- tmp$S_inv
  S     <- tmp$S
  
  
  # Initial values
  if ( is.null(init) ) {
    # Comment: Ideally make a function that generates initial values if the user does not provide any
    init <- generate_initial_values(n_region, nTimes, P, nKnots)#this function is added in utilities.R file
  }
  alpha        <- init$alpha
  beta         <- init$beta
  beta_tau     <- init$beta_tau
  eta          <- init$eta
  sigma2_alpha <- init$sigma2_alpha
  tau_i        <- init$tau
  tau0         <- init$tau0 
  w_i          <- init$w_i
  w0           <- init$w0


    
  # To store the MCMC samples
  nSave      <- nIter/nThin
  beta_mcmc  <- matrix( NA, dim(x)[2], nSave )
  alpha_mcmc <- array( NA, c(nTimes,n_region,nSave) )
  w_mcmc     <- array( NA, c(nKnots,n_region,nSave) )
  sigma_alpha_mcmc = rep(NA,nSave)
  tau_0_mcmc=rep(NA,nSave)
  beta_tau_mcmc=rep(NA,nSave)
  tau_i_mcmc = matrix(NA,length(init$tau),nSave)
  w_0_mcmc = matrix(NA,length(init$w0),nSave)
  
  
  
  # Comment: Remaining parameters could be added
  idx        <- 1 # Comment: this is to keep track where to store
  
  
  
  # Stepsize parameters 
  # Comment: **common** stepsize for all parameters is better, it is more efficient if they are updated all at once
  step             <- 0.005
  History          <- rep( NA, nBurnin ) 
  HistoryAccepted  <- rep( NA, nBurnin ) 
  HistoryAccepted0 <- rep( NA, nIter )     # Comment: this is for the entire chain, to track the acceptance rate in the end 
  
  
  
  # Batch choice
  tmp            <- minibatch_numbers( region, time, batch_type, prop, min_per_group )
  scaling_beta0  <- tmp$scaling_factor_beta
  scaling_w0     <- tmp$scaling_factor_w
  scaling_alpha0 <- tmp$scaling_factor_alpha
  indices        <- tmp$indices
  N_sample       <- tmp$N_sample
  N_total        <- tmp$N_total
  sample_size    <- sum(N_total)
  batch0         <- rep(1,length(region))
  
  
  
  # For the hybrid version
  batch_hybrid         <- batch0
  scaling_beta_hybrid  <- scaling_beta0*0 + 1
  scaling_w_hybrid     <- scaling_w0*0 +1
  scaling_alpha_hybrid <- scaling_alpha0*0 + 1
  
  
  
  ###########################################################################
  ###########################################################################
  # MCMC loop
  for ( iter in 1:nIter ) {
    
    # Select the batch
    if (batch_type!='full') {
      batch0 <- minibatch_select( indices, N_sample, N_total, sample_size )
    }
   
    
    # Adapt the batch and scaling factors for hybrid
    if ( hybrid & ((iter%%barker_frequency)==0) ) {
      batch    <- batch_hybrid
      sf_beta  <- scaling_beta_hybrid
      sf_w     <- scaling_w_hybrid
      sf_alpha <- scaling_alpha_hybrid
    } else {
      batch    <- batch0
      sf_beta  <- scaling_beta0
      sf_w     <- scaling_w0
      sf_alpha <- scaling_alpha0
    }
    
    
    # Update the random regression coefficients, random effects, spline parameters 
    #tmp                    <- barker_update( beta, alpha, w_i, step, x, dataset, Z, sigma2_alpha, S_inv, w0, tau_i, index,batch_index,scale_factor)
    tmp = barker_update_cpp( beta, alpha, w_i, step, Z, sigma2_alpha, S_inv, w0,  tau_i, x, index, y, batch, region, time, sf_beta, sf_w, sf_alpha )
    beta                   <- tmp$beta
    w_i                    <- tmp$w_i
    alpha                  <- tmp$alpha
    HistoryAccepted0[iter] <- 1*tmp$accepted
    
    
    #Gibbs step
    gibbs_out              <- gibbs_update( alpha_sigma, beta_sigma, alpha, w0, alpha_tau0, beta_tau0, S_inv, tau_i, shape_beta_tau, rate_beta_tau, shape_tau_i, w_i, S)
    sigma2_alpha           <- gibbs_out$ sigma2_alpha
    tau0                   <- gibbs_out$tau0
    beta_tau               <- gibbs_out$beta_tau
    tau_i                  <- gibbs_out$tau_i
    w0                     <- gibbs_out$w0
    
  
    # Adapt the stepsize 
    if ( iter<=nBurnin ) {
      History[iter]         <- step
      HistoryAccepted[iter] <- 1*tmp$accepted
      step                  <- adapt_stepsize( step, window, iter, target_acceptance, HistoryAccepted, History )
    }
    
    
    # If the burnin has just finished, make the stepsize the average of last burnin/2 iterations for the rest of the MCMC
    if ( iter==nBurnin ) {
      step <- mean( History[ (0.5*nBurnin + 1):nBurnin ] )
    }
    
    
    # Save the current state
    if ( (iter%%nThin)==0 ) {
      beta_mcmc[,idx]   <- beta
      alpha_mcmc[,,idx] <- alpha
      w_mcmc[,,idx]     <- w_i
      # Rest of parameters to be added here
      sigma_alpha_mcmc[idx] <-sigma2_alpha
      tau_0_mcmc[idx]<-tau0
      beta_tau_mcmc[idx]<-beta_tau
      tau_i_mcmc[,idx]<-tau_i
      w_0_mcmc[,idx]<-w0
      idx               <- idx + 1
      if ( iter%%500==0 ) {
        print(paste0('MCMC iteration ',iter))
      }
    }
  }
  ###########################################################################
  ###########################################################################
  
  
  
  
  # Return the MCMC samples 
  return( list(
    Z=Z,
    S=S,
    beta = beta_mcmc,
    alpha = alpha_mcmc,
    w = w_mcmc,
    sigma_alpha_mcmc=sigma_alpha_mcmc,
    tau_0_mcmc=tau_0_mcmc,
    beta_tau_mcmc=beta_tau_mcmc,
    tau_i_mcmc=tau_i_mcmc,
    w_0_mcmc=w_0_mcmc,
    index=index,
    step_history = History,
    acceptance_ratio = mean( HistoryAccepted0[-c(1:nBurnin)] ),  # only want to look at post-burnin acceptance
    acceptance = HistoryAccepted0, 
    N_sample = N_sample,
    sf_beta = scaling_beta0,
    sf_alpha = scaling_alpha0,
    sf_w = scaling_w0
  ))
}
  
