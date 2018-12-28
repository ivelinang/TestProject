#include "NormalDistribution1D.h"
#include <cmath>

//constructor
NormalDistribution1D::NormalDistribution1D(double mean, double stdev) :
	mean(mean), stdev(stdev)
{
	if (stdev <= 0.0)
	{
		throw("Standard deviation must be positive");
	}
	alpha = 1.0 / stdev;
	beta = -mean * alpha;
}

// destructor

NormalDistribution1D::~NormalDistribution1D()
{

}

// CDF(x)
double NormalDistribution1D::CDF(const double x) const
{
	return cody_cdf(alpha*x + beta);
}


// iCDF(p), i.e. CDF(x) = p
double NormalDistribution1D::iCDF(const double p) const
{
	return stdev * moro_icdf(p) + mean;
}

// PDF(x)
double NormalDistribution1D::PDF(const double x) const
{
	return pdf(alpha*x + beta);
}

double NormalDistribution1D::pdf(const double x)
{
	if (fabs(x) < 35.0) // to prevent exp() from underflowing
		return OneOverRootTwoPi * exp(-x * x / 2.0);
	else
		return 0.0;
}

// CDF Cordy W
// Cody W 1969 Rational Chebyshev approximations for the error function
// slow, very accurate, widely used

double NormalDistribution1D::cody_cdf(const double x)
{

	const double c[7] = { -0.0356098437018154,
		6.99638348861914,
		21.9792616182942,
		242.667955230532,
		15.0827976304078,
		91.1649054045149,
		215.058875869861 };

	const double xc[8] = { 0.000000136864857382717,
		     0.564195517478974,
		     7.21175825088309,
		    43.1622272220567,
		   152.989285046940,
		   339.320816734344,
		   451.918953711873,
		   300.459261020162 };

	const double yc[7] = { 12.7827273196294,
			77.0001529352295,
		   277.585444743988,
		   638.980264465631,
		   931.354094850610,
		   790.950925327898,
		   300.459260956983 };

	const double fx[5] = { 0.0223192459734185,
			 0.278661308609648,
			 0.226956593539687,
			 0.0494730910623251,
			 0.00299610707703542
	};

	const double gx[5] = { 1.98733201817135,
			 1.05167510706793,
			 0.191308926107830,
			 0.0106209230528468,
			 0.564189583547756 };

	double total, z, zz;

	z = (x >= 0) ? x : -x;
	z = z * OneOverRootTwo;

	if (z == 0.0)
		total = 0.5;
	else if (z < 0.46875)
	{
		zz = z * z;
		total = 0.5 - 0.5*z*(c[3] + zz * (c[2] + zz * (c[1] + zz * c[0]))) /
			(c[6] + zz * (c[5] + zz * (c[4] + zz)));
	}
	else if (z < 4.0)
	{
		zz = 1 / (z*z);
		total = 0.5*exp(-z * z)*(xc[7] + z * (xc[6] + z * (xc[5] + z * (xc[4] + z * (xc[3] + z * (xc[2] + z * (xc[1] + z * xc[0]))))))) /
			(yc[6] + z * (yc[5] + z * (yc[4] + z * (yc[3] + z * (yc[2] + z * (yc[1] + z * (yc[0] + z)))))));

	}
	else if (z < 25.0)
	{
		zz = 1 / (z*z);
		total = 0.5*exp(-z * z) / z * (OneOverRootTwoPi + zz * (-fx[4] + zz * (-fx[3] + zz * (-fx[2] + zz * (-fx[1] + zz * (-fx[0]))))) /
			(gx[3] + zz * (gx[2] + zz * (gx[1] + zz * (gx[0] + zz)))) 
			//+ gx[4]
			);
	}
	else
		total = 0.0;
	return (x >= 0) ? 1-total : total;

	


}

double NormalDistribution1D::moro_cdf(const double x)
{
	return 2;
}

double NormalDistribution1D::intfrac_cody_cdf(const double x, double *intpart, double *fracpart) {
	return 2;
}
double NormalDistribution1D::density_cody_cdf(const double x, double *density) {
	return 2;
}

double NormalDistribution1D::moro_icdf(const double p) {
	const double a[4] = {  2.50662823884,
		 -18.61500062529,
		 41.39119773534,
		 -25.44106049637 };

	const double b[4] = { -8.47351093090,
		23.08336743743,
		-21.06224101826,
		 3.13082909833 };

	const double c[9] = { 0.3374754822726147,
		 0.9761690190917186,
		 0.1607979714918209,
		 0.0276438810333863,
		 0.0038405729373609,
		 0.0003951896511919,
		 0.0000321767881768,
		 0.0000002888167364,
		 0.0000003960315187 };

	double total, t, tt, z;

	z = (p <= 0.5) ? p : 1 - p;

	if (z == 0.0)
	{
		total = -37.5;
	}
	else if (z > 0.08)
	{
		t = 0.5 - z;
		tt = t * t;
		total = -t * (a[0] + tt * (a[1] + tt * (a[2] + tt * (a[3])))) /
			(1.0 + tt * (b[0] + tt * (b[1] + tt * (b[2] + tt * (b[3] + tt * b[4])))));
	}
	else
	{
		t = log(-log(z));
		total = c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * (c[5] + t * (c[6] + t * (c[7] + t * c[8])))))));
		total = -total;
	}

	return (p <= 0.5) ? total : -total;
}
double NormalDistribution1D::hp_moro_icdf(const double x) {
	return 2;
}
double NormalDistribution1D::acklam_icdf(const double p) {
	return 2;
}