# Dissertation-Bayesian-MCMC


This project implements Bayesian inference using Markov Chain Monte Carlo (MCMC) methods with a focus on computational efficiency for large-scale data. \
Parameters with tractable full conditional distributions are updated using standard Gibbs sampling. The remaining parameters are updated using gradient-based MCMC methods.\
To further improve scalability, especially in the presence of large datasets, we replace standard gradient-based updates with Stochastic Gradient MCMC (SGMCMC) methods.
The workflow includes:

Defining a Bayesian hierarchical model\
Simulating data for validation\
Validating performance on simulated data\
Applying the model to real-world data for posterior inference and covariate effect estimation\
All code for model specification, algorithm implementation, and empirical evaluation is included in this repository.

