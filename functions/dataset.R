# Script that generates the dataset that we will use for the simulation

rm(list = ls())       


source('simulate_data.R')


# Dataset specifications 
prevalence  <- c(0.2,0.5,0.8)
sample_size <- c(30,90,300,900)
combination <- expand.grid( prevalence = prevalence, sample_size = sample_size)
nTimes      <- 10
sd_alpha    <- 0.1


# Generate the dataset
SEED <- 36531
set.seed(SEED)
dt   <- simulate_data( nTimes, combination$sample_size, sd_alpha, combination$prevalence )


# Save to use later
save(dt, file='simulated_dataset.RData')