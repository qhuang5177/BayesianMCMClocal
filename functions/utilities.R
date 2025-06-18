

# Function that sets up the spline design and penalisation matrix provided the number of time points and number of knots
spline_setup <- function( nTimes, nKnots )
{
  
  blank_data     <- data.frame(y = rnorm(nTimes, 0, 1), t = 1:nTimes)
  knots          <- list(t = seq(1, nTimes, length.out = nKnots))
  dummy_spline   <- jagam(y ~ s(t, bs = "cs", k = nKnots), file = "dummy.jags", data = blank_data)
  Z              <- dummy_spline$jags.data$X 
  S              <- dummy_spline$jags.data$S1
  S              <- rbind(0, cbind(0, S))
  S[1, 1]        <- 0.1                        # Precision for the intercept  
  Q              <- solve(S) 
  
  return(list( Z=Z, S_inv=S, S=Q ))#S:penalty matrix Q: inverse matrix
  
}


# For each individual provides of where to pick up eta_{it} and alpha_{it} from as.vector(Z%*%w) and as.vector(alpha), respectively 
# Comment: this only needs to be calculated once
individuals_index <- function( region, time )
{
  N        <- length(region)
  n_region <- length(table(region))
  nTimes   <- length(table(time))
  
  idx_region <- rep( 1:n_region, each=nTimes )
  idx_time   <- rep( 1:nTimes, times=n_region )
  n          <- nTimes*n_region
  
  index <- rep( NA, N )
  for ( i in 1:n ) {
    flag_region <- region==idx_region[i]
    flag_time   <- time==idx_time[i]
    index[ flag_region & flag_time ] = i
  }
  
  return(index)
}



# From vectorised Barker proposal, returns alpha, beta and w_i
from_vector <- function( params, nTimes, n_region, K, P )
{
  
  alpha_it_len = nTimes * n_region
  beta_len     = P
  w_i_len      = n_region * K
  
  # Random effects
  alpha = params[1:alpha_it_len]
  alpha = matrix(alpha, nrow = nTimes, ncol = n_region)
  
  # Regression coefficients
  beta = params[(alpha_it_len + 1):(alpha_it_len + beta_len)]
  
  # Spline parameters
  w_i = params[(alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]
  w_i = matrix(w_i, nrow = K, ncol = n_region)
  
  return(list( alpha=alpha, beta=beta, w_i=w_i))
}

#Flatten parameters
flatten_params <- function(alpha, beta, w_i) {
  return(c(as.vector(alpha), beta, as.vector(w_i)))
}







# Adapts the Barker stepsize
adapt_stepsize <- function( step, window, iter, target, HistoryAccepted, History )
{
  
  step_new <- step
  
  # Only do this if window iterations of the algorithm have run
  if ( iter>=window ) {
    
    # Look at the acceptance ratio in the past window iterations and adapt accordingly
    acceptance_ratio = mean( HistoryAccepted[(iter-window+1):iter] )
    if ( acceptance_ratio > target+0.05 ) {
      step_new <- step * 1.005
    } else if ( acceptance_ratio < target-0.05 ) {
      step_new <- step * 0.995
    }
    
  } 
  
  return( step_new )
  
}





# Initialize the parameters
generate_initial_values <- function(n_region, nTimes, P, nKnots) {
  return(list(
    alpha        = matrix(0, nrow = nTimes, ncol = n_region),       
    beta         = rep(0.1, P),                                    
    beta_tau     = 1,                                                
    eta          = matrix(0, nrow = nTimes, ncol = n_region),        
    sigma2_alpha = 1,                                                
    tau          = rep(2, n_region),                                 
    tau0         = 1,                                               
    w_i          = matrix(0.5, nrow = nKnots, ncol = n_region),      
    w0           = rep(0.5, nKnots)                                
  ))
}






