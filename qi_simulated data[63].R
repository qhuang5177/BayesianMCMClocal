nTimes = 14
N = 23
K = 6
S = diag(K)

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
