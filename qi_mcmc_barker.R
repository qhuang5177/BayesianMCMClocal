# Derivative of the log-posterior with respect to the parameters


nTimes=14
K=6
N=23
nPeople=100


gradient_alpha_it = function(alpha, sigma_alpha,y,p)
{
  alpha_it_gradient<-matrix(0,nTimes,N)
  for (i in 1:N) {
    for(t in 1:nTimes){
      grad_sum=0
      for (j in 1:nPeople) {
        grad_sum=grad_sum+(y[i,t,j]-p[i,t,j])
      }
     alpha_it_gradient[t,i]=grad_sum-alpha[t,i]/(sigma2_alpha^2)
    }
  }
  
return(as.vector(alpha_it_gradient))
  
}


gradient_beta=function(y,p,x,beta){
  
  beta_gradient=rep(0,K)
  grad_sum=0
  for (i in 1:N) {
    for(t in 1:nTimes){
      for (j in 1:nPeople) {
        grad_sum=grad_sum+(y[i,t,j]-p[i,t,j])*x[i,t,j]
      }
    
    }
    beta_gradient=grad_sum-(1/10)*beta
  }
  return(beta_gradient)
}

gradient_w_i=function(p,y,Z,S,tau_i,w_i,w_0){
  w_i_gradient=rep(0,K)
    for(t in 1:nTimes){
      for (j in 1:nPeople) {
        w_i_gradient=w_i_gradient+(y[i,t,j]-p[i,t,j])*Z[t,]
      }

    w_i_gradient=w_i_gradient-tau_i[i]*(solve(S)%*%(w_i[i]-w_0))
  }
  
  return(w_i_gradient)
  
}







