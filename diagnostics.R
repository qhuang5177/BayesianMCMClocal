

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













nSave = dim(fit$beta)[2]
nBurn = 1000
source('../plot.R')
p_vec = log(combination$prevalence/(1-combination$prevalence))

# Regression coefficients
P = dim(fit$beta)[1]
for ( i in 1:P ) {
  name = paste0( 'Covariate ', i) 
  plot(fit$beta[i,-c(1:nBurn)], main=name )
  abline( h=dt$beta[i] ,col=2, lwd=2, lty=2 )
}


# Smooth terms
smooth_mcmc <- array( NA, c(nSave,nTimes,I) )
for ( i in 1:nSave ) {
  smooth_mcmc[i,,] = fit$Z %*% fit$w[,,i]
}
for ( i in 1:I ) {
  name = paste0( 'Smooth term region ', i) 
  polygons( smooth_mcmc[-c(1:nBurn),,i], MAIN=name)
  abline( h=p_vec[i] ,col=2, lwd=2, lty=2, main=name )
}


# Random effects
alpha_mcmc <- array( NA, c(nSave,nTimes,I) )
for ( i in 1:nSave ) {
  alpha_mcmc[i,,] = fit$alpha[,,i]
}
for ( i in 1:I ) {
  name = paste0( 'Random effects region ', i) 
  boxes( alpha_mcmc[-c(1:nBurn),,i], MAIN=name, outliers=FALSE )
  points( dt$alpha[,i], pch=19, col=2, cex=2 )
}



# Combined 
p_mcmc = alpha_mcmc+smooth_mcmc
for ( i in 1:I ) {
  name = paste0( 'Logit-prevalence region ', i) 
  polygons( p_mcmc[-c(1:nBurn),,i], MAIN=name)
  lines( p_vec[i]+dt$alpha[,i], pch=19, col=2, cex=2 )
}
