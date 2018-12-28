#pragma once

#include <boost/noncopyable.hpp>

/*
calculate pdf and cdf of bivariate normal distribution
 mean mu = [m1, mu2] 
 variances = [sigma1^2, sigma2^2]
 correlation matrix = corr = [[1, rho[, [rho,1]]

 Three methods:
  - Drezner
  - Drezner & Wesolowsky
  - Divgi

  Drezner Z. "Computation of the Bivariate Normal Integral",
  Mathematics of computation 32 1878

  Drezner Z Wesolowsky G "On the computation of the Bivariate
  Normal Integral", Journal of Statist. Comp. Simul 1989

  Divgi D.R "Calculation of univariate and bivariate normal probability functions"
  Annals of Statistics 1979

  For comparison b/n three methods see

  Genz A. "Numerical computation of rectangular bivariate and trivariate
  normal and t probabilities"

  SPEED and ACCURACY

  NAME			ACCURACY (in dp)			RELATIVE SPEED
  ---------------------------------------------------------
  Drezner			5							0.5
  DreznerWeslowsky
	plain			7							1.0
  Divgi
	plain			8							1.5
	hp				11							1.2

 dp - absolute error is bounded by 1e-7

 Drezner: not very reliable, especially for rho close to 1

 Drezner_wesolowsky: very solid for values of rho close to 1, 
 but has a ditch in the region around a=-0.3, b=-0.3, rho = -0.6

 Divgi: Very solid for values of rho close to 1. 
 Divgi performs slightly better regarding speed and accuracy

*/

class NormalDistribution2D : boost::noncopyable
{
public:
	explicit NormalDistribution2D(const double rho);
	explicit NormalDistribution2D(const double mean1 = 0.0, const double mean2 = 0.0,
		const double stdev1 = 1.0, const double stdev2 = 1.0,
		const double rho = 0.0);
	~NormalDistribution2D();

	double PDF(const double x, const double y) const;
	double CDF(const double x, const double y) const;

	static double pdf(const double x, const double y, const double rho);
	static double divgi_cdf(const double x, const double y, const double rho);

private:
	double alpha1, beta1;
	double alpha2, beta2;
	const double rho;

};