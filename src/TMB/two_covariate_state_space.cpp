#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_INTEGER(nyears);
  DATA_VECTOR(y);
  DATA_MATRIX(Xij);
  DATA_MATRIX(Uij);
  DATA_IVECTOR(yearindex);
  DATA_VECTOR(x1tobs);
  DATA_IVECTOR(x1yearindex);
  DATA_SCALAR(nx1);
  DATA_VECTOR(x2tobs);
  DATA_IVECTOR(x2yearindex);
  DATA_SCALAR(nx2);
  DATA_SCALAR(priorlogsigmaobs1_mu);
  DATA_SCALAR(priorlogsigmaobs1_sigma);
  DATA_SCALAR(priorlogsigmaobs2_mu);
  DATA_SCALAR(priorlogsigmaobs2_sigma);
  
  PARAMETER_VECTOR(beta); // fixed effect coefficients
  PARAMETER_VECTOR(gamma); // random effect coefficients
  PARAMETER(logsigma_gamma); // standard deviation of random effects
  PARAMETER_VECTOR(vt); // the vector of vt's
  PARAMETER(logitrho);
  PARAMETER(logsigma_vt); // variation in r_t
  PARAMETER(logno); // initial log population size
  PARAMETER(logphi); // negative binomial dispersion parameters
  PARAMETER_VECTOR(vx1t); // vector of the vxt's
  PARAMETER(logsigma_vx1t);
  PARAMETER(logsigma_x1obs);
  PARAMETER(x1o);
  PARAMETER(thetax1);
  
  PARAMETER_VECTOR(vx2t); // vector of the vxt's
  PARAMETER(logsigma_vx2t);
  PARAMETER(logsigma_x2obs);
  PARAMETER(x2o);
  PARAMETER(thetax2);
  
  Type jnll = 0.0;
  jnll -= dnorm(logsigma_x1obs, priorlogsigmaobs1_mu,priorlogsigmaobs1_sigma, true);
  jnll -= dnorm(logsigma_x2obs, priorlogsigmaobs2_mu,priorlogsigmaobs2_sigma, true);
  
  
  
  // Transformed parameters
  
  Type sigma_gamma = exp(logsigma_gamma);
  Type sigma_v = exp(logsigma_vt);
  Type phi = exp(logphi);
  Type rho = 1 / (1 + exp( logitrho ) );
  Type nyears_minus_1 = nyears - 1;
  vector<Type> rt(nyears - 1);
  vector<Type> lognt(nyears);
  vector<Type> nt(nyears);
  
  // Transformed parameters for covariate
  Type sigma_vx1t = exp(logsigma_vx1t);
  Type sigma_x1obs = exp(logsigma_x1obs);
  vector<Type> x1t(nyears);
  
  Type sigma_vx2t = exp(logsigma_vx2t);
  Type sigma_x2obs = exp(logsigma_x2obs);
  vector<Type> x2t(nyears);
  
  ADREPORT(sigma_gamma);
  ADREPORT(sigma_v);
  ADREPORT(phi);
  ADREPORT(rho);
  ADREPORT(sigma_vx1t);
  ADREPORT(sigma_vx2t);
  ADREPORT(sigma_x1obs);
  ADREPORT(sigma_x2obs);
  
  
  // State space model for x1
  
  x1t(0) = x1o;
  
  for ( int i  = 1; i < nyears; i++){ 
    x1t(i) = x1t(i - 1) + vx1t(i - 1);
  }
  
  for (int i = 0; i < nx1; i++) {
    jnll -= dnorm(x1tobs(i), x1t(x1yearindex(i)), sigma_x1obs,true);
  }
  ADREPORT(x1t);
  
  
  
  // repeat for covariate 2
  
  x2t(0) = x2o;
  
  for ( int i  = 1; i < nyears; i++){ 
    x2t(i) = x2t(i - 1) + vx2t(i - 1);
  }
  
  for (int i = 0; i < nx2; i++) {
    jnll -= dnorm(x2tobs(i), x2t(x2yearindex(i)), sigma_x2obs,true);
  }
  ADREPORT(x2t);
  
  
  rt(0) = vt(0) + thetax1 * x1t(0) + thetax2 * x2t(0);
  lognt(0) = logno;
  nt(0) = exp(lognt(0));
  
  for ( int i  = 1; i < nyears_minus_1; i++){ 
    rt(i) = rt(i - 1) * rho + pow(1 - pow(rho, 2), -2) * vt(i) + thetax1 * x1t(i) + thetax2 * x2t(i);
  }
  
  for ( int i  = 1; i < nyears; i++){ 
    lognt(i) = lognt(i - 1) + rt(i - 1);
    nt(i) = exp(lognt(i));
  }
  ADREPORT(nt);
  ADREPORT(rt);
  ADREPORT(lognt);
  
  
  // generate a vector of lognt that align with observations
  vector<Type> logntdata(n);
  for (int i =0; i < n; i++) {
    logntdata(i) = lognt(yearindex(i));
  }
  
  
  // random effects likelihood
  jnll -= sum(dnorm(gamma,0,sigma_gamma, true));
  jnll -= sum(dnorm(vt, 0, sigma_v, true)); 
  jnll -= sum(dnorm(vx1t, 0, sigma_vx1t, true)); 
  jnll -= sum(dnorm(vx2t, 0, sigma_vx2t, true)); 
  
  vector<Type>  mu(n);
  mu = exp( logntdata + Xij* beta + Uij * gamma) ;
  
  
  vector<Type> var(n);
  for(int i = 0; i < n; i++){ 
    var[i] = mu[i] + pow(mu[i],2) * pow (phi, -1);
    jnll -=  	dnbinom2(y[i], mu[i],var[i], true);
    
  }
  
  return jnll;
}
