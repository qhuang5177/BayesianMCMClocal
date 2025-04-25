#load('Y.RData')
#load('X.RData')

#Y <- sim_data$y

#X <- as.matrix(sim_data[, c("x1", "x2", "x3")])

set.seed(123)
library(MASS)

# set true value
#true_sigma2_alpha <- 0.1
#true_tau_0 <- 0.5
#true_beta_tau <- 1
#true_shape_tau_i <- 2.1
#true_rate_tau_i <- true_beta_tau


#alpha_it_new <- matrix(rnorm(14 * 23, mean = 0, sd = sqrt(true_sigma2_alpha)), nrow = 14, ncol = 23)
#tau_i_new <- rgamma(23, shape = true_shape_tau_i, rate = true_rate_tau_i)
#tau_0_new <- true_tau_0
#w_0_new <- mvrnorm(1, mu = rep(0, 6), Sigma = solve((1/tau_0_new) * S))




load('initial_values.RData')

S<-diag(K)


source('qi_mcmc_main[49].R')
posterior = stoch_gradient_logistic (nIters=12000, nThin=10, initial_values, alpha_sigma=10, beta_sigma=1,alpha_tau0 = 0.1,beta_tau0 = 0.1,shape_beta_tau=1,rate_beta_tau=1,shape_tau_i=2.1,tau_i_rate=1.3)
plot(1:1200,posterior$sigma_alpha_mcmc,type = 'l',main = 'Trace Plot of Sigma Random Effect ',xlab='Iterations',ylab='')
abline(h=0.1,col='red')
plot(1:1200,posterior$tau_0_mcmc,type = 'l',main = expression("Trace Plot of " * tau[0]),xlab='Iterations',ylab = '')
abline(h=0.7,col='red')

plot(1:1200,posterior$beta_tau_mcmc,type = 'l',main = expression("Trace Plot of " * beta[tau]),xlab='Iterations',ylab = '')
abline(h=1.3,col='red')

par(mfrow=c(4,4))
for (i in 1:16) {
  plot(1:1200,posterior$tau_i_mcmc[,i],type = 'l',main = expression("Trace Plot of " * tau[i]),xlab='Iterations',ylab = '')
  abline(h=tau_i[i],col='red',lwd=2)
}

par(mfrow=c(4,2))
for (i in 17:23) {
  plot(1:1200,posterior$tau_i_mcmc[,i],type = 'l',main = expression("Trace Plot of " * tau[i]),xlab='Iterations',ylab = '')
  abline(h=tau_i[i],col='red',lwd=2)
}


par(mfrow=c(3,2))
for (i in 1:6) {
  plot(1:1200,posterior$w_0_mcmc[,i],type = 'l',main = expression("Trace Plot of " * omega[i]),xlab='Iterations',ylab = '')
  abline(h=w_0[i],col='red',lwd=2)
}





