# 加载必要的包
library(MASS)
library(mgcv)
library(ggplot2)



set.seed(42)


nTimes <- 14    
nKnots <- 6   
blank_data <- data.frame(y = rnorm(nTimes, 0, 1), t = 1:nTimes)

knots <- list(t = seq(1, nTimes, length.out = nKnots))
dummy_spline <- jagam(y ~ s(t, bs = "cs", k = nKnots), file = "dummy.jags", data = blank_data)
Z <- dummy_spline$jags.data$X 


S <- dummy_spline$jags.data$S1
S <- rbind(0, cbind(0, S))
S[1, 1] <- 0.1
S <- solve(S)   

### 2. 模型总体设置
I <- 23         # Number of region
T <- nTimes     # Number of year
N <- 10000   # Number of individuals in ith region and jth year 
p <- 3          # Number of covariates

### simulate hyperparameter
beta_tau <- rgamma(1, shape = 1, rate = 1)        # Gamma(1,1)
tau0 <- rgamma(1, shape = 0.1, rate = 0.1)          # Gamma(0.1,0.1)
# w0 ~ Normal(0, (1/tau0)*S)，
w0 <- mvrnorm(1, mu = rep(0, nKnots), Sigma = (1/tau0)*S)

# beta ~ Normal(0,10I) 
beta <- mvrnorm(1, mu = rep(0, p), Sigma = diag(10, p))

# sigma_alpha^2 ~ InverseGamma(0.1,0.1)

sigma2_alpha <- 1/rgamma(1, shape = 10, scale = 1)
sigma_alpha <- sqrt(sigma2_alpha)

### Simulate f_it = eta_it + alpha_it for every region

# dimension :23*14
f_mat <- matrix(NA, nrow = I, ncol = T)

#do it vectorized and get T calculate next level 
for(i in 1:I){
  # For every region i，模拟 tau_i ~ Gamma(2.1, beta_tau)
  tau_i <- rgamma(1, shape = 5, rate = beta_tau)
  # simulate  w_i ~ Normal(w0, (1/tau_i)*S)
  w_i <- mvrnorm(1, mu = w0, Sigma = (1/tau_i)*S)
  eta_i <- as.vector(Z %*% w_i)
  # For yea t, simulate alpha_it ~ Normal(0, sigma2_alpha)
  alpha_i <- rnorm(T, mean = 0, sd = sigma_alpha)
  f_mat[i, ] <- eta_i + alpha_i
}


sim_data <- data.frame(
  region = rep(1:I, each = T * N),
  time = rep(rep(1:T, each = N), times = I)
)


# x1 ~ Bernoulli(0.5), x2 ~ Uniform(-1,1), x3 ~ Normal(0,1)
sim_data$x1 <- rbinom(n = nrow(sim_data), size = 1, prob = 0.5)
sim_data$x2 <- runif(n = nrow(sim_data), min = -1, max = 1)
sim_data$x3 <- rnorm(n = nrow(sim_data), mean = 0, sd = 1)

# Distribute the random effect term to every individual
# 对于每个区域 i 和年份 t，f_mat[i,t] 对应于该组合下的截距
sim_data$f <- NA
for(i in 1:I){
  for(t in 1:T){
    idx <- which(sim_data$region == i & sim_data$time == t)
    sim_data$f[idx] <- f_mat[i, t]
  }
}

### 6. 根据模型生成每个个体的响应

sim_data$lp <- sim_data$f + beta[1]*sim_data$x1 + beta[2]*sim_data$x2 + beta[3]*sim_data$x3

#  expit function
expit <- function(x) { 1 / (1 + exp(-x)) }
sim_data$p <- expit(sim_data$lp)

# Simulate y ~ Bernoulli(p)
sim_data$y <- rbinom(n = nrow(sim_data), size = 1, prob = sim_data$p)

T <- 14

par(mar=c(5, 4, 4, 4) + 0.3)  # 增加右边边距
plot(1:T, eta_i, type="l", col="blue", ylab="eta_i", xlab="Time (Year)",
     main="eta_i (left axis) and alpha_i (right axis) changing the mean of tao")
par(new=TRUE)
plot(1:T, alpha_i, type="l", col="red", axes=FALSE, xlab="", ylab="")
axis(side=4)
mtext("alpha_i", side=4, line=3)
legend("topleft", legend=c("eta_i", "alpha_i"), col=c("blue", "red"), lty=1)

#f_it <- eta_i + alpha_i
#plot(1:T, f_it, type="l", col="green", main="f_it = eta + alpha", ylab="f_it", xlab="Year",lwd=2)
#lines(1:T,eta_i,co='blue',lwd=2,lty=2)





hist(sim_data$p, breaks = 50, col = "lightblue", main = "Prior Predictive of p change the mean of tao", xlab = "p")

table(sim_data$y) / nrow(sim_data)  # y = 0 / 1 的比例


