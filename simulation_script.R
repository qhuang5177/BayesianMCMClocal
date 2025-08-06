

# generate the dataset
rm(list = ls())     


# libraries and scripts
source('./functions/simulate_data.R')
library(mgcv)
source('./functions/utilities.R')
source('./functions/mcmc.R')
source('./functions/gibbs.R')
library(Rcpp)
sourceCpp('./functions/barker.cpp')


# Dataset specifications 
prevalence  <- c(0.2,0.5,0.8)
sample_size <- c(30,90,300,900)
combination <- expand.grid( prevalence = prevalence, sample_size = sample_size)
nTimes      <- 10
sd_alpha    <- 0.1


# Generate the dataset
ID <- 1
SEED <- 5177+ID-1
set.seed(SEED)
dt   <- simulate_data( nTimes, combination$sample_size, sd_alpha, combination$prevalence )
dataset <- dt$sim_data
I       <- length(table(dataset$region))


# Simulation specifications
nThin   <- 200
nBurnin <- nThin*200
nMCMC   <- nThin*1500
nIter   <- nMCMC+nBurnin
nKnots  <- 3


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


# Full MCMC 
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="full"  
  ) 
})
results1 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)


# Stratified
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="stratified", prop=0.1, MH_correction=0, step_correction=0.5    
  ) 
})
results2 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)


# Stratified, different stepsize
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="stratified", prop=0.1, MH_correction=0, step_correction=0.75    
  ) 
})
results3 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)


# Stratified 20%
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="stratified", prop=0.2, MH_correction=0, step_correction=0.5    
  ) 
})
results4 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)


# Stratified with minimum sample size 
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="stratified_min", prop=0.1, MH_correction=0, step_correction=0.5, min_per_group=30    
  ) 
})
print(fit$N_sample)
results5 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)


# Stratified with minimum sample size 
mcmc_time <- system.time({
  fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data, covariates=c(3:8), nThin=nThin,
                      target_acceptance=0.55, batch_type="stratified", prop=0.1, MH_correction=0, step_correction=0.5, hybrid=T, barker_frequency=10
  ) 
})
results6 <- mcmc_results( fit, mcmc_time, nIter, nBurnin/nThin, dt)



# Save all the outputs
save( results1, file=paste0('./posteriors/res1_',ID,'.RData') )
save( results2, file=paste0('./posteriors/res2_',ID,'.RData') )
save( results3, file=paste0('./posteriors/res3_',ID,'.RData') )
save( results4, file=paste0('./posteriors/res4_',ID,'.RData') )
save( results5, file=paste0('./posteriors/res5_',ID,'.RData') )
save( results6, file=paste0('./posteriors/res6_',ID,'.RData') )





