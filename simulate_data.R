# Libraries
library(MASS)
library(mgcv)
library(ggplot2)


# Function that generates datasets for simulation studies/testing
# Comment: Best practice is to make this a function rather than call the script
# Comment: Made sample size to be a vector rather than a single number
simulate_data = function(
    nTimes,              # number of time points
    nKnots,              # number of knots
    N,                   # vector of sample size(reigion*how many people withiin it)
    p_vec = NULL
    )
{
  
  
  # Spline design matrices
  blank_data     <- data.frame(y = rnorm(nTimes, 0, 1), t = 1:nTimes)
  knots          <- list(t = seq(1, nTimes, length.out = nKnots))
  dummy_spline   <- jagam(y ~ s(t, bs = "cs", k = nKnots), file = "dummy.jags", data = blank_data)
  Z              <- dummy_spline$jags.data$X 
  S              <- dummy_spline$jags.data$S1
  S              <- rbind(0, cbind(0, S))
  S[1, 1]        <- 0.1                        # Precision for the intercept  
  S              <- solve(S)   
  
  
  # Simulate hyperparameters
  #beta_tau <- rgamma(1, shape =1, rate = 100)                        # Gamma(1,1)
  #tau0     <- rgamma(1, shape =10, rate =1)                          # Gamma(0.1,0.1)
  #w0       <- mvrnorm(1, mu = rep(0, nKnots), Sigma = ((1/tau0)*S))  # Prior for mean spline weights
  
  
  # Regression coefficients: beta ~ Normal(0, 10I)
  p    <- 6 # Comment: defined the number of variables here
  # Comment: variance of the regression coefficients too large. It makes sense to have a weekly informative prior, but generating from it creates extreme values of prevalence
  beta <- mvrnorm(1, mu = rep(0, p), Sigma = (diag(.01, p)))
 

  # Variance of alpha_it 
  #sigma2_alpha <- 1/rgamma(1, shape =10, rate = 1)
  #sigma_alpha  <- sqrt(sigma2_alpha)
  
  
  # Region-specific parameters
  I <- length(N)
  # Comment: better to save all these quantities so you can check if they have been estimated okay
  f_mat <- eta <- alpha <- matrix( NA, nrow = I, ncol=nTimes )
  p_vec <- rep(p_vec, length.out = I)
  
  #w     <- matrix( NA, nrow=I, ncol=nKnots ) 
  #tau   <- rep(NA,I)
  #for(i in 1:I){
    #tau[i]     <- rgamma(1, shape = 10, rate = beta_tau)                 # tau_i ~ Gamma
   # w[i,]      <- mvrnorm(1, mu = w0, Sigma = ( (1/tau[i])*S ) )         # w_i ~ Normal(w0, precision)
    #eta[i,]    <- as.vector( Z %*% w[i,] )                               # η_i = Z * w_i
    #alpha[i,]  <- rnorm( nTimes, mean = 0, sd =0.1)             # α_it ~ N(0, σ²)
    #f_mat[i, ] <- eta[i,] + alpha[i,]                                    # f_it = η_it + α_it
  #}
  
  
  for(i in 1:I){
    eta[i, ]   <- log(p_vec[i]/(1-p_vec[i]))
    alpha[i, ] <- rnorm(nTimes, mean = 0, sd = 0.08)
    f_mat[i, ] <- eta[i, ] + alpha[i, ]
  }
  

  # Prepare simulated dataset
  # Comment: this has changed a bit because the sample size changes over regions. Could extend to matrix N with sample size changing over time and units
  region <- time <- NULL
  for ( i in 1:I ) {
    region <- c( region, rep(i,nTimes*N[i]) )
    time   <- c( time, rep(1:nTimes,each=N[i]) )
  }
  sim_data        <- data.frame( region=region, time=time )
  
  
  # Simulate covariates
  # Comment: would be great if this could be more generall, e.g. user specifies how many covariates and/or types
  # Comment: better to define a quantity that is used multiple times, here nrow(sim_data)
  # Comment: perhaps one reason that prevalences are extreme is that variance of covariates is high?
  n           <- nrow(sim_data)
  sim_data$x1 <- rbinom(n, size = 1, prob = 0.5)     # x1 ~ Bernoulli
  sim_data$x2 <- runif(n, min = -1, max = 1)         # x2 ~ Uniform
  sim_data$x3 <- rnorm(n, mean = 0, sd = 1)    
  sim_data$x4 <- rnorm(n, mean = 0, sd = 2)    
  sim_data$x5 <- rnorm(n, mean = 1, sd = 1)    
  sim_data$x6 <- rnorm(n, mean = 2, sd = 2.4)
  
  
  # Assign random intercept f_it to each individual
  # Comment: I'd say avoid using 'T' for names of variables
  sim_data$f <- NA
  for(i in 1:I){
    idx1 <- sim_data$region == i
    # Comment: in the original code you are calculating this more times than necessary. When the sample size is large, it recude speed. This is a generic comment, not so much for this specific task
    for(t in 1:nTimes){
      idx2                   <- sim_data$time == t
      sim_data$f[idx1&idx2]  <- f_mat[i, t]
    }
  }
  
  
  # Generate linear predictor and probabilities
  # Comment: would be ideal to have more general code for the line below, e.g. put all X into a matrix?
  # exclude_names <- c("region", "time", "f", "lp", "p", "y")
  #covariate_names <- setdiff(names(sim_data), exclude_names)
  #X_mat <- as.matrix(sim_data[, covariate_names])
  #sim_data$lp <- sim_data$f + as.vector(X_mat %*% beta)
  
  # Comment: there was a typo in the previous code, all last 3 terms were beta[4]*x4
  sim_data$lp <- sim_data$f + beta[1]*sim_data$x1 + beta[2]*sim_data$x2 + beta[3]*sim_data$x3+beta[4]*sim_data$x4+beta[5]*sim_data$x5+beta[6]*sim_data$x6
 
   expit       <- function(x) { 1 / (1 + exp(-x)) }
  sim_data$p  <- expit(sim_data$lp)
  
  
  # Simulate y ~ Bernoulli(p)
  sim_data$y <- rbinom(n = nrow(sim_data), size = 1, prob = sim_data$p)

  
  # Return the dataset to the user
  # Comment: good to return a list to the user with simulated values of all parameters so you can test the code
  # Comment: in this code w is KxI whereas in the MCMC it is IxK
  return(
    list(
      sim_data     = sim_data, 
      beta         = beta,
      eta          = t(eta),
      alpha        = t(alpha),
      Z            = Z
    
    )
  )

}















