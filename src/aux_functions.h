#include <RcppArmadillo.h>
using namespace Rcpp;
double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double truncn2(double mu, double sigma, double lower, double upper);
arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma);