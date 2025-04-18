load('Y.RData')
load('X.RData')

Y <- sim_data$y

X <- as.matrix(sim_data[, c("x1", "x2", "x3")])

set.seed(34324)
#testing temporarily 
library(MASS)

# set true value
true_sigma2_alpha <- 0.1
true_tau_0 <- 0.5
true_beta_tau <- 1
true_shape_tau_i <- 2.1
true_rate_tau_i <- true_beta_tau


alpha_it_new <- matrix(rnorm(14 * 23, mean = 0, sd = sqrt(true_sigma2_alpha)), nrow = 14, ncol = 23)
tau_i_new <- rgamma(23, shape = true_shape_tau_i, rate = true_rate_tau_i)
tau_0_new <- true_tau_0
w_0_new <- mvrnorm(1, mu = rep(0, 6), Sigma = solve((1/tau_0_new) * S))





load('initial_values.RData')


S = diag(1,6)
source('qi_mcmc_main.R')
posterior = stoch_gradient_logistic (Y,X,nIters=12000, nThin=10, initial_values, alpha_sigma=10, beta_sigma=1,alpha_tau0 = 0.1,beta_tau0 = 0.1,shape_beta_tau=1,rate_beta_tau=1,shape_tau_i=2.1,tau_i_rate=beta_tau_new)
plot(1:1200,posterior$tau_i_mcmc,type = 'l')

mean(posterior$sigma_alpha_mcmc)
mean(posterior$tau_0_mcmc)
mean(posterior$beta_tau_mcmc)
mean(posterior$tau_i_mcmc)



