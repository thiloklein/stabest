#include <cmath> // used for INFINITY and isinf
#include "truncnorm_approx.h"
//#include "myrandom.h"

double qnorm_approx(double p)
{
	// see Abramowitz and Stegun, 26.2.23, http://people.math.sfu.ca/~cbm/aands/page_933.htm
	// |approx error| <= 4.5*10^-4
	double t, x_p;
	if (p>=1) return INFINITY;
	else if (p<=0) return -INFINITY;
	else {
	t = p<0.5 ? sqrt(log(1/(p*p))) : sqrt(log(1/((1-p)*(1-p))));
	x_p =  t - ((0.010328*t +0.802853)*t + 2.515517) / (((0.001308*t + 0.189269)*t + 1.432788)*t + 1.0);
	return p<0.5 ? -x_p : x_p;
	}
}

double pnorm_approx(double x)
{
	// see Abramowitz and Stegun, 26.2.17, http://people.math.sfu.ca/~cbm/aands/page_932.htm
	// |approx error| <= 7.5*10^-8
	double p_x, t;
	if (std::isinf(x))
		return x>0 ? 1.0 : 0.0;
	t = 1 / (1 + 0.2316419*fabs(x));
	p_x = exp(-x*x/2) / sqrt(2*3.141593) * t * (0.31938153 + t*(-0.356563782 + t*(1.781477937 + t*(-1.821255978 + t*1.330274429 ))));
	return x>=0 ? 1-p_x : p_x;
}

// Uniform random number generator, see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
/*
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	
double runif()
{
	std::uniform_real_distribution<double> dis(0, 1);
	return dis(gen);
}
*/

double truncnorm_approx(double p, double mean, double sigma, double lower, double upper)
{
  // draws from the inverse cdf of a truncated normal variable, i.e. with
  // p \in [0,1] we it returns x ~ N(0,1) s.th. lower <= x <= upper and Pr(X<=x) = p
	double u,x;
	// transform lower and upper bounds to N(0,1)
	lower = pnorm_approx((lower - mean) / sigma);
	upper = pnorm_approx((upper - mean) / sigma);
	u = lower + p * (upper - lower);
	x = qnorm_approx(u);
	return x*sigma + mean;
}



