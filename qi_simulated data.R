# Load required libraries
library(MASS)
library(mgcv)
library(ggplot2)

set.seed(42)

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

# 2. Overall model settings
I <- 23         # Number of regions
T <- nTimes     # Number of years
N <- 10 # Number of individuals per region-year
p <- 6   # Number of covariates

# Simulate hyperparameters
beta_tau <- rgamma(1, shape =1, rate = 100)            # Gamma(1,1)
tau0 <- rgamma(1, shape =10, rate =1)            # Gamma(0.1,0.1)
w0 <- mvrnorm(1, mu = rep(0, nKnots), Sigma = ((1/tau0)*S))  # Prior for mean spline weights

# Regression coefficients: beta ~ Normal(0, 10I)
beta <- mvrnorm(1, mu = rep(0, p), Sigma = (diag(10, p)))

# Variance of alpha_it 
sigma2_alpha <- 1/rgamma(1, shape =10, rate = 1)
sigma_alpha <- sqrt(sigma2_alpha)

# Simulate f_it = eta_it + alpha_it for each region
f_mat <- matrix(NA, nrow = I, ncol = T)

for(i in 1:I){
  tau_i <- rgamma(1, shape = 10, rate = beta_tau)                 # tau_i ~ Gamma
  w_i <- mvrnorm(1, mu = w0, Sigma = ((1/tau_i)*S))                # w_i ~ Normal(w0, precision)
  eta_i <- as.vector(Z %*% w_i)                                  # η_i = Z * w_i
  alpha_i <- rnorm(T, mean = 0, sd = sigma_alpha)                # α_it ~ N(0, σ²)
  f_mat[i, ] <- eta_i + alpha_i                                  # f_it = η_it + α_it
}

# Prepare simulated dataset
sim_data <- data.frame(
  region = rep(1:I, each = T * N),
  time = rep(rep(1:T, each = N), times = I)
)

# Simulate covariates
sim_data$x1 <- rbinom(n = nrow(sim_data), size = 1, prob = 0.5)     # x1 ~ Bernoulli
sim_data$x2 <- runif(n = nrow(sim_data), min = -1, max = 1)         # x2 ~ Uniform
sim_data$x3 <- rnorm(n = nrow(sim_data), mean = 0, sd = 1)    
sim_data$x4 <- rnorm(n = nrow(sim_data), mean = 0, sd = 2)    
sim_data$x5 <- rnorm(n = nrow(sim_data), mean = 1, sd = 1)    
sim_data$x6 <- rnorm(n = nrow(sim_data), mean = 2, sd = 2.4)    

# x3 ~ Normal

# Assign random intercept f_it to each individual
sim_data$f <- NA
for(i in 1:I){
  for(t in 1:T){
    idx <- which(sim_data$region == i & sim_data$time == t)
    sim_data$f[idx] <- f_mat[i, t]
  }
}

# Generate linear predictor and probabilities
sim_data$lp <- sim_data$f + beta[1]*sim_data$x1 + beta[2]*sim_data$x2 + beta[3]*sim_data$x3+beta[4]*sim_data$x4+beta[4]*sim_data$x4+beta[4]*sim_data$x4
expit <- function(x) { 1 / (1 + exp(-x)) }
sim_data$p <- expit(sim_data$lp)

# Simulate y ~ Bernoulli(p)
sim_data$y <- rbinom(n = nrow(sim_data), size = 1, prob = sim_data$p)

# Visualization
par(mar=c(5, 4, 4, 4) + 0.3)
plot(1:T, eta_i, type="l", col="blue", ylab="eta_i", xlab="Time (Year)",
     main="eta_i (left axis) and alpha_i (right axis)")
par(new=TRUE)
plot(1:T, alpha_i, type="l", col="red", axes=FALSE, xlab="", ylab="")
axis(side=4)
mtext("alpha_i", side=4, line=3)
legend("topleft", legend=c("eta_i", "alpha_i"), col=c("blue", "red"), lty=1)

# Histogram of prior predictive p
hist(sim_data$p, breaks = 50, col = "lightblue",probability = T, 
     main = "Prior Predictive of p", xlab = "p")


