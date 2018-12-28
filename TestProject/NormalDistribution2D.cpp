

#include "NormalDistribution1D.h"
#include "NormalDistribution2D.h"
#include "Constants.h"

NormalDistribution2D::NormalDistribution2D(const double rho) : alpha1(1.0), alpha2(1.0),
beta1(0.0), beta2(0.0), rho(rho)
{
	if (rho*rho >= 1)
	{
		throw("Correlation is not in the interval [-1,1]");
	}
}

//general constructor

NormalDistribution2D::NormalDistribution2D(const double mean1,
	const double mean2, const double stdev1, const double stdev2,
	const double rho) :rho(rho)
{
	if (rho*rho >= 1)
	{
		throw("Correlation is not in the interval [-1,1]");
	}

	if (stdev1*stdev2 <= 0)
	{
		throw("standard deviations must be positive ");
	}

	alpha1 = 1.0 / stdev1;
	beta1 = -mean1 * alpha1;

	alpha2 = 1.0 / stdev2;
	beta2 = -mean2 * alpha2;
}

//Destructor

NormalDistribution2D::~NormalDistribution2D()
{}

// pdf
double NormalDistribution2D::PDF(const double x, const double y) const
{
	return pdf(alpha1*x + beta1, alpha2*y + beta2, rho);
}

double NormalDistribution2D::CDF(const double x, const double y) const
{
	return divgi_cdf(alpha1*x + beta1, alpha2*y + beta2, rho);
}

double NormalDistribution2D::pdf(const double a, const double b, const double rho)
{
	double inner, invdet;
	invdet = 1 / (1 - rho * rho);
	inner = 0.5*((a - rho * b)*a + (b - rho * a)*b)*invdet;

	if (inner < 25.0)
		return OneOverTwoPi * sqrt(invdet)*exp(-inner);
	else
		return 0.0;
}

double NormalDistribution2D::divgi_cdf(double a, double b, const double rho)
{
	const int N = 11; //chosen as Divgi paper
	const double TOLERANCE = 1e-15; 

	static const double d[N + 1] = {
		1.253309, -0.999907, 0.6260210,
		-0.3310819, 0.1518857, -0.6004582e-1,
		0.1977542e-1, -0.5151399e-2, 0.9963604e-3,
		-0.132076766e-3, 0.105772e-4, -0.3826e-6 };

	return 2;
}