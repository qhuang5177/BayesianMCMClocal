nTimes = 14
N = 23
K = 6
S = diag(K)


library(MASS)
library(mgcv)
library(ggplot2)



nTimes <- 14    
nKnots <- 6   

# Dummy dataset for generating spline basis
blank_data <- data.frame(y = rnorm(nTimes, 0, 1), t = 1:nTimes)
knots <- list(t = seq(1, nTimes, length.out = nKnots))
dummy_spline <- jagam(y ~ s(t, bs = "cs", k = nKnots), file = "dummy.jags", data = blank_data)
Z <- dummy_spline$jags.data$X 

S <- dummy_spline$jags.data$S1
S <- rbind(0, cbind(0, S))
S[1, 1] <- 0.1  
S <- solve(S)   

S<-diag(6)





sigma_alpha = 0.1
alpha = matrix( rnorm(N*nTimes,0,sqrt(sigma_alpha)),nTimes,N)
tau_0 = 0.7
w_0 = mvrnorm( 1, rep(0,K), S/tau_0 )
beta_tau = 1.3
tau_i = rgamma( N, shape=2.1, rate=beta_tau )
w_i = matrix(NA, K, N)
for ( i in 1:N ) {
  w_i[,i] = mvrnorm(1,w_0,S/tau_i[i])
}

eta_it <- matrix(NA, nrow = nTimes, ncol = N)
for (i in 1:N) {
  eta_it[, i] <- Z %*% w_i[, i]
}

f_it <- eta_it + alpha  


for (i in 1:N) {
  for (t in 1:nTimes) {
    for (j in 1:J) {
      x_vec <- X_itj[i, t, j, ]  # length K
      linear_pred <- f_it[t, i] + sum(beta * x_vec)
      p <- 1 / (1 + exp(-linear_pred))  # expit
      p_itj[i, t, j] <- p
      y_itj[i, t, j] <- rbinom(1, 1, p)
    }
  }
}
















initial_values = list(
  sigma_alpha=0.1,
  alpha = alpha,
  tau_0 = tau_0,
  w_0=w_0,
  beta_tau=beta_tau,
  tau_i=tau_i,
  w_i=w_i
)

save(initial_values,file='initial_values.RData')
