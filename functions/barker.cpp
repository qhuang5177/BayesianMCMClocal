


// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace arma;



/* Select batch */ 
// [[Rcpp::export]]
arma::uvec minibatch_select_cpp(const arma::umat &idx_range, const arma::uvec &idx_vec, unsigned int N, const arma::uvec &N_sample, const arma::uvec &N_total)
{
  
  unsigned int C = idx_range.n_cols;
  arma::uvec batch(N, fill::zeros);
  
  for (unsigned int c=0; c<C; c++) {
    
    unsigned int n = N_sample(c);
    if (n == 0) continue;  
    
    unsigned int idx0          = idx_range(0, c);
    unsigned int idx1          = idx_range(1, c);
    unsigned int total_indices = N_total(c);
    
    /* Fill the indices to sample from */
    arma::uvec indices(total_indices);
    for (unsigned int k=0; k<total_indices; k++) {
      indices(k) = idx_vec(idx0 + k);
    }
    
    // Fisher-Yates shuffle for sampling without replacement
    for (unsigned int i = 0; i < n; i++) {
      unsigned int j = i + std::floor(R::runif(0.0, 1.0) * (total_indices - i));
      std::swap(indices(i), indices(j));
      batch(indices(i)) = 1;
    }
  }
  
  return batch;
}



/* Put alpha, beta and w in a single vector */ 
// [[Rcpp::export]]
vec to_vector( const vec &beta, const mat &alpha, const mat &w )
{
  unsigned int P = beta.n_elem;
  unsigned int N = alpha.n_cols;
  unsigned int T = alpha.n_rows;
  unsigned int K = w.n_rows;
  unsigned int i, j, idx = 0;
  
  vec parameter(P+N*T+K*N);
  
  /* Random effects */ 
  for ( i=0 ; i<N ; i++ ) {
    for ( j=0 ; j<T ; j++ ) {
      parameter(idx) = alpha(j,i);
      idx += 1;
    }
  }
  
  /* Regression coefficients */ 
  for ( i=0 ; i<P ; i++ ) {
    parameter(idx) = beta(i);
    idx += 1;
  }
  
  /* Spline coefficients */ 
  for ( i=0 ; i<N ; i++ ) {
    for ( j=0 ; j<K ; j++ ) {
      parameter(idx) = w(j,i);
      idx += 1;
    }
  }
  
  return(parameter);
}



/* Save the proposed values */ 
void from_vector( const vec &proposed, vec &beta, mat &alpha, mat &w )
{
  unsigned int D = proposed.n_elem;
  unsigned int N = alpha.n_cols;
  unsigned int T = alpha.n_rows;
  unsigned int K = w.n_rows;
  unsigned int P = beta.n_elem;
  unsigned int i, j, idx=0;
  
  /* Random effects */
  for ( i=0 ; i<N ; i++ ) {
    for ( j=0 ; j<T ; j++ ) {
      alpha(j,i) = proposed(idx);
      idx += 1;
    }
  }
  /* Regression coefficients */
  for ( i=0 ; i<P ; i++ ) {
    beta(i) = proposed(idx);
    idx += 1;
  }
  /* Spline coefficients */
  for ( i=0 ; i<N ; i++ ) {
    for ( j=0 ; j<K ; j++ ) {
      w(j,i) = proposed(idx);
      idx += 1;
    }
  }
  
}



/* Function to evaluate the prevalence */
vec p_fitted( const vec &beta, const mat &w, const mat &alpha, const mat &x, const mat &Z, const uvec &index, vec &prob, const uvec &batch, unsigned int MH )
{
  mat eta = Z * w;
  vec f   = vectorise(eta+alpha);
  unsigned int N = x.n_rows;
  vec linear_predictor(N,fill::zeros);
  unsigned int i;
  if ( MH==1 ) {
    /* When there is MH correction */
    /* We evaluate linear predictor for everyone, as needed  for the likelihood */ 
    /* Fitted probabilities only for those in the batch */ 
    for ( i=0 ; i<N ; i++ ) {
      linear_predictor(i) = dot(x.row(i).t(),beta) + f(index(i)-1);
      if ( batch(i)==1 ) {
        prob(i) = 1.0/( 1.0+exp(-linear_predictor(i)));
      }
    }
  } else {
    /* When there is no MH correction */
    /* Only need fitted probabilities (and hence linear predictors) for those in the batch */
    for ( i=0 ; i<N ; i++ ) {
      if ( batch(i)==1 ) {
        linear_predictor(i) = dot(x.row(i).t(),beta) + f(index(i)-1);
        prob(i) = 1.0/( 1.0+exp(-linear_predictor(i)));
      }
    }
  }
  
  /* Return the linear predictor */ 
  return(linear_predictor);
}



/* Log-posterior evaluation */ 
double log_post(const uvec &y, const vec &lp, const vec &beta, double beta_var, const mat &alpha, 
                double sigma_alpha, const mat &w, const vec &w0, const mat &S_inv, const vec &tau )
{
  double c;
  
  /* Likelihood contributions */
  c = sum( y % lp - log(1.0 + exp(lp)) );
  
  /* Prior contributions */ 
  c += -0.5*dot(beta,beta)/beta_var;
  c += -0.5*accu(alpha%alpha)/sigma_alpha;
  unsigned int n = w.n_cols;
  unsigned int K = w.n_rows;
  vec tmp(K);
  for ( unsigned int i=0 ; i<n ; i++ ) {
    tmp = w.col(i)-w0;
    c += -0.5*tau(i)*dot( tmp, S_inv*tmp );
  }
  
  
  return(c);
}



/* Function to evaluate gradients */
vec grad( const uvec &y, const vec &p, const mat &x, const vec &beta, double beta_var, const uvec &batch,
            const uvec &region, const uvec &time, const mat &alpha, double sigma_alpha, const mat &Z, const mat &w,
            const vec&tau, const mat &S_inv, const vec &w0, double scaling_beta, const vec &scaling_w,
            const mat &scaling_alpha )
{
  
  /* Regression coefficients */ 
  double tmp;
  unsigned int i, j,  N = y.n_elem;
  unsigned int P = beta.n_elem;
  vec grad_beta(P,fill::zeros);
  /* likelihood */
  for ( i=0 ; i<N ; i++ ) {
    if ( batch(i)==1 ){
      tmp = -p(i) + y(i);
      for ( j=0 ; j<P ; j++ ) {
        grad_beta(j) += tmp*x(i,j);
      }
    }
  }
  /* prior and scaling factor */
  for ( j=0 ; j<P ; j++ ) {
    grad_beta(j) /= scaling_beta;
    grad_beta(j) -= beta(j)/beta_var;
  }
  
  
  /* Random effects */ 
  unsigned int n = alpha.n_cols;
  unsigned int T = alpha.n_rows;
  mat grad_alpha(T,n,fill::zeros);
  /* Likelihood */
  for ( i=0 ; i<N ; i++ ) {
    if ( batch(i)==1 ){
      tmp = -p(i) + y(i);
      grad_alpha( time(i)-1, region(i)-1 ) += tmp;  
    }
  }
  /* scaling factor and prior */
  for ( i=0 ; i<n ; i++ ) {
    for ( j=0 ; j<T ; j++ ) {
      if ( scaling_alpha(j,i)>0 ) {
        grad_alpha(j,i) /= scaling_alpha(j,i);
      }
      grad_alpha( j, i ) -= alpha(j,i)/sigma_alpha;
    }
  }
  
  
  /* Spline coefficients */ 
  unsigned int K = w.n_rows; 
  mat grad_w(K,n,fill::zeros);
  /* likelihood */
  for ( i=0 ; i<N ; i++ ) {
    if ( batch(i)==1 ){
      tmp = -p(i) + y(i);
      for ( j=0 ; j<K ; j++ ) {
        grad_w( j, region(i)-1 ) += tmp*Z(time(i)-1,j);
      }
    }
  }
  /* scaling factor and prior */ 
  for ( i=0 ; i<n ; i++){
    if ( scaling_w(i)>0 ) {
      grad_w.col(i) /= scaling_w(i);
    }
    grad_w.col(i) -= tau(i)*( S_inv*(w.col(i)-w0) );
  }
  
    
  return( join_cols( vectorise(grad_alpha), grad_beta, vectorise(grad_w) ) );
}




// [[Rcpp::export]]
Rcpp::List barker_update_cpp( const vec &beta, const mat &alpha, const mat &w, double stepsize, const mat &Z, 
                              double sigma_alpha, const mat &S_inv, const vec &w0, const vec &tau, const mat &x,
                              const uvec &index, const uvec &y, const uvec &batch, 
                              const uvec &region, const uvec &time, double scaling_beta, const vec &scaling_w,
                              const mat &scaling_alpha, unsigned int MH
)
{
  
  /* Current value */
  vec current = to_vector( beta, alpha, w);
  
  
  /* Noise */
  unsigned int D = current.n_elem;
  double noise_sd = 0.1*stepsize;
  vec noise(D);
  unsigned int i;
  for ( i=0 ; i<D ; i++ ) {
    noise(i) = R::rnorm( stepsize, noise_sd );
  }
  
  
  /* Fitted probabilities current value */ 
  unsigned int n = x.n_rows;
  vec p_current(n);
  vec lp_current = p_fitted(beta, w, alpha, x, Z, index, p_current, batch, MH);
  
  
  /* Evaluate the gradient at the current point */
  vec beta_x = grad( y, p_current, x, beta, 10.0, batch, region, time, alpha, sigma_alpha, Z, w, tau, S_inv, w0, scaling_beta, scaling_w, scaling_alpha);
  
  
  /* Propose a new value */ 
  vec b(D,fill::zeros);
  vec q(D,fill::zeros);
  double u;
  for ( i=0 ; i<D ; i++ ) {
    q(i) = 1.0/(1.0+exp(-beta_x(i)*noise(i)));
    u = R::runif(0.0,1.0);
    if ( u<=q(i) ) {
      b(i) = 1.0;
    } else {
      b(i) = -1.0;
    }
  }
  vec proposed = current + b % noise;
  
  
  /* Extract the proposed values */
  unsigned int P = beta.n_elem;
  unsigned int N = alpha.n_cols;
  unsigned int T = alpha.n_rows;
  unsigned int K = w.n_rows;
  vec beta_new(P);
  mat alpha_new(T,N);
  mat w_new(K,N);
  from_vector( proposed, beta_new, alpha_new, w_new );
  
  
  /* Gradient and fitted probabilities proposed value */ 
  vec p_proposed(n);
  vec lp_proposed = p_fitted(beta_new, w_new, alpha_new, x, Z, index, p_proposed, batch, MH );
  vec beta_y = grad( y, p_proposed, x, beta_new, 10.0, batch, region, time, alpha_new, sigma_alpha, Z, w_new, tau, S_inv, w0, scaling_beta, scaling_w, scaling_alpha);
  
  
  /* MH correction */ 
  double R;
  if ( MH==1 ) {
    /* Likelihoods */
    double post_current  = log_post(y, lp_current, beta, 10.0, alpha, sigma_alpha, w, w0, S_inv, tau );
    double post_proposed = log_post(y, lp_proposed, beta_new, 10.0, alpha_new, sigma_alpha, w_new, w0, S_inv, tau );
    R = post_proposed-post_current;
    /* Correction */
    vec tmp1  = -beta_y % (current-proposed);
    vec tmp2  = -beta_x % (proposed-current);
    vec zero_vec = arma::zeros<vec>(tmp1.n_elem);
    vec term1 = arma::max(tmp1, zero_vec) + arma::log1p(arma::exp(-arma::abs(tmp1)));
    vec term2 = arma::max(tmp2, zero_vec) + arma::log1p(arma::exp(-arma::abs(tmp2)));
    double correction = arma::sum(-term1 + term2);
    R += correction;
  } else {
    R = 1.0;
  }
  

  /* Accept-reject */ 
  int accepted = 0;
  u = log(R::runif(0.0,1.0));
  if ( u<=R ) {
    accepted = 1;
  } else {
    beta_new = beta;
    alpha_new = alpha;
    w_new = w;
  }

  
  
  /* Output */ 
  return Rcpp::List::create(
    Rcpp::Named("grad") = beta_x,
    Rcpp::Named("log_post") = R,
    //Rcpp::Named("b") = b,
    //Rcpp::Named("noise") = noise,
    //Rcpp::Named("current") = current
    Rcpp::Named("beta") = beta_new,
    Rcpp::Named("alpha") = alpha_new,
    Rcpp::Named("w_i") = w_new,
    Rcpp::Named("accepted") = accepted
  );
}