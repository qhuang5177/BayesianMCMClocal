minibatch_numbers = function( region, time, method, prop=0.2, min_sample=1 )
{
  
  # For faster sampling of indices in SG, we will save the ID by region and time
  idx     = 1:length(region)
  indices = list()
  
  # First calculate the sample size per time and region
  nRegions <- length(table(region))
  nTimes   <- length(table(time)) 
  N        <- matrix(NA,nTimes,nRegions) 
  for ( i in 1:nRegions ) {
    indices[[i]] = list()
    idx1 = region==i
    for ( t in 1:nTimes ) {
      idx2 = time==t
      indices[[i]][[t]] = idx[idx1&idx2]
      N[t,i] = sum(idx1&idx2)
    }
    names(indices)
  }
  
  
  # Now calculate how many individuals to sample per region
  if ( method=='full' ) {
    N_sample = N
  } else if ( method=='stratified' ) {
    N_sample = round(N*prop)
    N_sample[N_sample<=0] = 1
  } else if ( method=='stratified_min' ) {
    N_sample = 0*N+min_sample
    N_sample = (N<N_sample)*N + (N>=N_sample)*N_sample
    N_rest = round(length(region)*prop)-sum(N_sample)
    if (N_rest<=0) {
      print(paste0('With the specified minimum, the proportion of data used is ',round(100*sum(N_sample)/length(region),2),'%'))
    } else {
      N_sample = N_sample + round(N*(N_rest/length(region)))
    }
  } else {
    stop("The possible methods are 'full/stratified/stratified_min'")
  }
  
  
  # Now find the scaling factors
  scaling_factor_beta  = sum(N_sample)/length(region)
  scaling_factor_w     = apply(N_sample,2,sum)/apply(N,2,sum)
  scaling_factor_alpha = N_sample/N 
  
  
  # output
  return(list(
    scaling_factor_beta=scaling_factor_beta,
    scaling_factor_w=scaling_factor_w,
    scaling_factor_alpha=scaling_factor_alpha,
    N_total=N,
    N_sample=N_sample,
    indices=indices
  ))
  
}



minibatch_select = function( indices, N_sample, N_total, sample_size )
{
  nRegions = dim(N_total)[2]
  nTimes   = dim(N_total)[1]
  batch = rep(0,sample_size)
  for (i in 1:nRegions) {
    for (j in 1:nTimes) {
      batch[sample(indices[[i]][[j]],N_sample[j,i],replace=FALSE)]=1
    }
  }
  return(batch)
}


# Function that sets up the spline design and penalisation matrix provided the number of time points and number of knots
spline_setup <- function( nTimes, nKnots )
{
  
  blank_data     <- data.frame(y = rnorm(nTimes, 0, 1), t = 1:nTimes)
  knots          <- list(t = seq(1, nTimes, length.out = nKnots))
  dummy_spline   <- jagam(y ~ s(t, bs = "cs", k = nKnots), file = "dummy.jags", data = blank_data)
  Z              <- dummy_spline$jags.data$X 
  S              <- dummy_spline$jags.data$S1
  S              <- rbind(0, cbind(0, S))
  S[1, 1]        <- 0.1                        # Precision for the intercept  
  Q              <- solve(S) 
 
  return(list( Z=Z, S_inv=S, S=Q ))#S:penalty matrix Q: inverse matrix
  
  
}


# For each individual provides of where to pick up eta_{it} and alpha_{it} from as.vector(Z%*%w) and as.vector(alpha), respectively 
# Comment: this only needs to be calculated once
individuals_index <- function( region, time )
{
  N        <- length(region)
  n_region <- length(table(region))
  nTimes   <- length(table(time))
  
  idx_region <- rep( 1:n_region, each=nTimes )
  idx_time   <- rep( 1:nTimes, times=n_region )
  n          <- nTimes*n_region
  
  index <- rep( NA, N )
  for ( i in 1:n ) {
    flag_region <- region==idx_region[i]
    flag_time   <- time==idx_time[i]
    index[ flag_region & flag_time ] = i
  }
  
  return(index)
}



# From vectorised Barker proposal, returns alpha, beta and w_i
from_vector <- function( params, nTimes, n_region, K, P )
{
  
  alpha_it_len = nTimes * n_region
  beta_len     = P
  w_i_len      = n_region * K
  
  # Random effects
  alpha = params[1:alpha_it_len]
  alpha = matrix(alpha, nrow = nTimes, ncol = n_region)
  
  # Regression coefficients
  beta = params[(alpha_it_len + 1):(alpha_it_len + beta_len)]
  
  # Spline parameters
  w_i = params[(alpha_it_len + beta_len + 1):(alpha_it_len + beta_len + w_i_len)]
  w_i = matrix(w_i, nrow = K, ncol = n_region)
  
  return(list( alpha=alpha, beta=beta, w_i=w_i))
}

#Flatten parameters
flatten_params <- function(alpha, beta, w_i) {
  return(c(as.vector(alpha), beta, as.vector(w_i)))
}







# Adapts the Barker stepsize
adapt_stepsize <- function( step, window, iter, target, HistoryAccepted, History )
{
  
  step_new <- step
  
  # Only do this if window iterations of the algorithm have run
  if ( iter>=window ) {
    
    # Look at the acceptance ratio in the past window iterations and adapt accordingly
    acceptance_ratio = mean( HistoryAccepted[(iter-window+1):iter] )
    if ( acceptance_ratio > target+0.05 ) {
      step_new <- step * 1.005
    } else if ( acceptance_ratio < target-0.05 ) {
      step_new <- step * 0.995
    }
    
  } 
  
  return( step_new )
  
}





# Initialize the parameters
generate_initial_values <- function(n_region, nTimes, P, nKnots) {
  return(list(
    alpha        = matrix(0, nrow = nTimes, ncol = n_region),       
    beta         = rep(0, P),                                    
    beta_tau     = 1,                                                
    eta          = matrix(0, nrow = nTimes, ncol = n_region),        
    sigma2_alpha = 1,                                                
    tau          = rep(2, n_region),                                 
    tau0         = 1,                                               
    w_i          = matrix(0, nrow = nKnots, ncol = n_region),      
    w0           = rep(0.5, nKnots)                                
  ))
}


#batch sampler
select_batch <- function(dataset, 
                         batch_type = c("random", "stratified"), 
                         prop = 0.2, 
                         min_per_group = 10) {
  
  # Step 1: 计算总目标样本数
  N <- nrow(dataset)
  batch_size <- ceiling(N * prop)
  
  # 创建 region-time 分组 index
  dataset$group <- interaction(dataset$region, dataset$time, drop = TRUE)
  
  # Step 2: 每个 group 至少抽 min_per_group
  grouped_data <- split(dataset, dataset$group)
  selected_list_1 <- lapply(grouped_data, function(group_data) {
    group_n <- nrow(group_data)
    n_sample <- min(group_n, min_per_group)
    group_data[sample(1:group_n, n_sample,replace = F), ]
  })
  sample_data_1 <- do.call(rbind, selected_list_1)
  
  # Step 3: 从剩余样本中补足 remaining_n
  remaining_n <- batch_size - nrow(sample_data_1)
  if (remaining_n <= 0) {
    rownames(sample_data_1) <- NULL
    return(sample_data_1)
  } else {
    # 找出未被抽中的 observation
    used_ids <- rownames(sample_data_1)
    remaining_data <- dataset[!rownames(dataset) %in% used_ids, ]
    
    if (batch_type == "random") {
      sample_data_2 <- remaining_data[sample(1:nrow(remaining_data), remaining_n,replace = F), ]
      
    } else if (batch_type == "stratified") {
      # stratify by region
      grouped_remain <- split(remaining_data, remaining_data$region)
      group_sizes <- sapply(grouped_remain, nrow)
      group_prop <- group_sizes / sum(group_sizes)
      group_n_sample <- round(group_prop * remaining_n)
      
      sample_data_2 <- do.call(rbind, lapply(names(grouped_remain), function(g) {
        dat <- grouped_remain[[g]]
        n <- min(nrow(dat), group_n_sample[g])
        dat[sample(1:nrow(dat), n,replace = F), ]
      }))
    }
    
    final_sample <- rbind(sample_data_1, sample_data_2)
    rownames(final_sample) <- NULL
    return(final_sample)
  }
} 








