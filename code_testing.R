# Clear the workspace
rm(list = ls())       # Comment: always best to start a simulation script with this so that no objects from previous sessions corrupt any calculations



# libraries and scripts
library(mgcv)
#source('./function0718/simulate_data.R')
source('/Users/qhuang2024/Library/Mobile Documents/com~apple~CloudDocs/Dissertation cam  /BayesianMCMClocal/function0718/mcmc.R')
source('/Users/qhuang2024/Library/Mobile Documents/com~apple~CloudDocs/Dissertation cam  /BayesianMCMClocal/function0718/simulate_data.R')
source('/Users/qhuang2024/Library/Mobile Documents/com~apple~CloudDocs/Dissertation cam  /BayesianMCMClocal/function0718/barker.R')
source('/Users/qhuang2024/Library/Mobile Documents/com~apple~CloudDocs/Dissertation cam  /BayesianMCMClocal/function0718/gibbs.R')
source('/Users/qhuang2024/Library/Mobile Documents/com~apple~CloudDocs/Dissertation cam  /BayesianMCMClocal/function0718/utilities.R')




# Simulate a dataset
nTimes <- 5
#I      <- 6
#N      <- rep(100,I)


nKnots <- 3
prevalence = c(0.2,0.5,0.8)
sample_size = c(20,50,100)

combination <- expand.grid(prevalence = prevalence,
                sample_size = sample_size)


I<-n_region <- nrow(combination)  
N <- combination$sample_size



dt=simulate_data( nTimes=nTimes, nKnots=nKnots, N=N,p_vec=prevalence)
#dt=simulate_data( nTimes=nTimes, nKnots=nKnots, N=N)






# Set the MCMC specifications and initial values
nBurnin <-300000
nIter   <- 2*nBurnin
nThin   <- 100 
# Comment: could make a function to generate initial values. Using the values you used to simulate the data makes testing easy
#init    <- list(
 # alpha        = 0*dt$alpha,
  #beta         = 0*dt$beta,
  #beta_tau     = dt$beta_tau,
  #eta          = dt$eta,
  #sigma2_alpha = dt$sigma2_alpha,
  #tau          = dt$tau,
  #tau0         = dt$tau0,
  #w_i          = 0*dt$w,
  #w0           = dt$w0
#)

#I modified the code for simulate the data. So now I didnot use previous way to generate initial value
init <- list(
  alpha        = 0 * dt$alpha,
  beta         = 0 * dt$beta,
  beta_tau     = 0.01,                 
  eta          = dt$eta,            
  sigma2_alpha = 0.05,                
  tau          = rep(1, I),         
  tau0         = 10,
  w_i          = t(matrix(0, I, nKnots)),
  w0           = rep(0, nKnots)
)








# Run the MCMC
mcmc_time <- system.time({
fit <- barker_mcmc( nIter=nIter, nBurnin=nBurnin, init=init, nKnots=nKnots, dataset=dt$sim_data,target_acceptance=0.55, covariates=c(3:8), nThin=nThin, alpha_sigma = 10,        #default value for hyperparpmeter
                    beta_sigma = 1,
                    alpha_tau0 = 10, 
                    beta_tau0 =  1,
                    shape_beta_tau = 1, 
                    rate_beta_tau = 100,
                    shape_tau_i = 10,             # Thinning: if there are too many MCMC iterations, you can't save all of the samples
                    use_batch =TRUE, 
                    mix_full_barker=TRUE,
                    barker_frequency = 20, 
                    batch_type = "random",   # Default to random batch if using batch
                    prop = 0.2,              # Proportion of sample to be drawn if using batch
                    min_per_group = 20,
                    window=100              # Minimum samples per group for stratified sampling
) 
})

# 输出耗时（单位：秒）
cat("Total time taken: ", round(mcmc_time["elapsed"], 2), "seconds\n")








# Look at the stepsize
plot(fit$step_history,type='l',xlab='MCMC iteration',ylab='Stepsize')
print(fit$acceptance_ratio) # acceptance rate


  # Find the posterior of the prevalence
nSave           <- dim(fit$beta)[2]
prevalence_mcmc <- array( NA, c(nSave,nTimes,I) )



for (i in 1:nSave){
  prevalence_mcmc[i,,] <- fit$Z%*%fit$w[,,i] + fit$alpha[,,i]
}
expit       <- function(x) { 1 / (1 + exp(-x)) }



#prevalence_true <- expit(dt$Z%*%dt$w + dt$alpha)
# For constant prevalence，obtained by using C$prevalence
prevalence_true <- expit(dt$eta + dt$alpha)   # (nTimes × I)

# Make some plots/bias compare with the true value
test_region=3
test_time=5
plot(expit(prevalence_mcmc[,test_time,test_region]),type='l',ylim=c(0,1),xlab='MCMC iteration',ylab='prevalence')
abline(h=prevalence_true[test_time,test_region],col=2,lwd=2,lty=2)







library(coda)

#calculate ESS
ess_matrix <- matrix(NA, nrow = nTimes, ncol = I)

for (t in 1:nTimes) {
  for (r in 1:I) {
    samples <- expit(prevalence_mcmc[, t, r])
    ess_matrix[t, r] <- effectiveSize(samples)
  }
}


colnames(ess_matrix) <- paste0("Region_", 1:I)
rownames(ess_matrix) <- paste0("Time_", 1:nTimes)


print(round(ess_matrix, 2))










#calculating time
mcmc_sec <- as.numeric(mcmc_time["elapsed"])
mcmc_min <- round(mcmc_sec / 60, 2)


#
 save(fit,file='fit_SGMCMC_mix_600000.RData')


 
 
 
 
 
 
 
 
 
 
 
 
 
