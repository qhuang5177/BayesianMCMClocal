# Clear the workspace
rm(list = ls())       # Comment: always best to start a simulation script with this so that no objects from previous sessions corrupt any calculations
ID = 1


# libraries and scripts
library(mgcv)
source('utilities.R')
source('mcmc.R')
source('utilities.R')
source('gibbs.R')
library(Rcpp)
sourceCpp('barker.cpp')



# Load the dataset
load('simulated_dataset.RData')
dataset <- dt$sim_data
I       <- length(table(dataset$region))



# Specifications of knots and MCMC iterations
nKnots  <- 3
nThin   <- 250
nBurnin <- nThin*500
nMCMC   <- nThin*2000
nIter   <- nMCMC+nBurnin



# Initial values
init <- list(
  alpha        = 0 * dt$alpha,
  beta         = 0 * dt$beta,
  beta_tau     = 0.01,                 
  eta          = 0 * dt$eta,            
  sigma2_alpha = 0.05,                
  tau          = rep(1, I),         
  tau0         = 10,
  w_i          = t(matrix(0, I, nKnots)),
  w0           = rep(0, nKnots)
)



# Run the MCMC
mcmc_time <- system.time({
fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data,target_acceptance=0.55, covariates=c(3:8), nThin=nThin, 
                    alpha_sigma      = 0.01,        
                    beta_sigma       = 0.01,
                    alpha_tau0       = 10, 
                    beta_tau0        = 1,
                    shape_beta_tau   = 1, 
                    rate_beta_tau    = 100,
                    shape_tau_i      = 10,            
                    window           = 100, 
                    batch_type       = "stratified",     
                    hybrid           = FALSE,          
                    barker_frequency = 10,    # How often to do Barker in hybrid
                    prop = 0.1,               # Proportion of sample to be drawn if using batch
                    min_per_group = 10        # Minimum samples per group for stratified sampling
                    ) 
})
name <- paste0('./posteriors/regular_',ID,'.RData')
save(fit,file=name)

# 输出耗时（单位：秒）
cat("Total time taken: ", round(mcmc_time["elapsed"], 2), "seconds\n")




# Save the MCMC fit
save(fit,file='reference_chain.RData')





