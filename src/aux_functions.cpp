// Some helper functions - se stabit2Sel2.cpp
#define ARMA_NO_DEBUG
//#define DEBUG
#include <RcppArmadillo.h>
#include "aux_functions.h"
using namespace Rcpp;


// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

double norm_rs(double a, double b)
{
  double  x;
  x = Rf_rnorm(0.0, 1.0);
  while( (x < a) || (x > b) ) x = norm_rand();
  return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

double half_norm_rs(double a, double b)
{
  double   x;
  x = fabs(norm_rand());
  while( (x<a) || (x>b) ) x = fabs(norm_rand());
  return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// ======================================================================

double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) xstar = 0.0;
  else xstar = a;
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

double exp_rs(double a, double b)
{
  double  z, u, rate;
  
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a)) z = R::rexp(rate);
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a)) z = R::rexp(rate);
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

// truncn2( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

double truncn2(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf) || (b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  // Second scenario
  else if((a * b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2) z = unif_rs(a,b);
    else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
    else z = exp_rs(a,b);
    if(change) z = -z;
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

// ---------------------------------------------  
// random multivariate normal sample generator using RcppArmadillo
// from http://gallery.rcpp.org/articles/simulate-multivariate-normal/
// ---------------------------------------------  

arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma) {
#ifdef DEBUG
  printf("mvrnormArma: mu is %d x %d, \n", mu.n_rows, mu.n_cols);
#endif
  arma::rowvec y = arma::randn(1,mu.n_rows); // as<arma::rowvec>(arma::randn(1,ncols)); //arma::randn(1,ncols)
  // include error handling for cholesky inverse - sigma must be positive definite,
  // see https://math.stackexchange.com/questions/462682/why-does-the-cholesky-decomposition-requires-a-positive-definite-matrix
  arma::mat sigma_chol = arma::zeros(sigma.n_rows, sigma.n_rows); // empty container for output (square)
  bool success = false;
  int Ntries = 100; // number of tries, be generous on this
  while (success==false && Ntries>=0) {
    success = arma::chol(sigma_chol, sigma); // does not raise a runtime exception
    if (!success) {
      if (!Ntries) throw std::runtime_error( "chol(): decomposition failed" ); //
      Rcpp::Rcout << "cholesky inversion failed, with matrix sigma:" << std::endl;
      Rcpp::Rcout << sigma << std::endl;
      Rcpp::Rcout << "adding 1e-6 to diagonal elements:" << std::endl;
      sigma += arma::eye(sigma.n_rows,sigma.n_rows) * 1e-6;
      Ntries -= 1;
    }
  }
  return arma::trans( arma::trans(mu) + y*sigma_chol);
}


