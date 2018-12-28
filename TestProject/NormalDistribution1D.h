#pragma once

#include <boost/noncopyable.hpp>

/*
The function in this class computes N(x), the CDF of a normally distributed random variable with mean = mu
and variance = sigma^2

Function to compute density n(x) is also included

Four algorightms:
1. An old formula of unknown source: old_normal_cdf
2. C. Hastings formula
3. B. Moros
4. W. Cody

Possible qualifiers:
density_   computes PDF and CDF at the same time
intfrac_   the result is split in integral and fractional part
asm_       implementation is in assembly code available (faster)
sse_       implementation in assembly with SSE instructions available (much faster)
pair_      computes CDF() of two values at the same time


SPEED and ACCURACY:

NAME				ACCURACY				RELATIVE SPEED
-------------------------------------------------------------
Cody
	plain				15						1.0
	density				15						1.0
	intfrac				15						1.0
	density+intfrac		15						0.9

Moro
	plain				9						2.4
	intfrac				9						2.4
	asm					9						3.0

Hastings
	"slow"
	plain				8						1.2
	density				8						1.2
	intfrac				8						1.2
	density+intfrac		8						1.2

	"fast"
	plain				7						1.9	
	intfrac				7						1.9
	sse					6						2.4
	sse+pair			6						4.9

Old
	plain				7						0.9

Accuracy: 15dp menas absolute error is bounded by 1e-15
Speed ratio compared to Cody algorithm
Speed is tested using Intel Pentium III



*/

class NormalDistribution1D : boost::noncopyable
	{
	public:
		NormalDistribution1D(double mean = 0.0, double variane = 1.0);
		~NormalDistribution1D();

		// computed PDF(x) of N(mean, stdev)
		double PDF(const double x) const;
		// computed CDF(x) of N(mean, stdev)
		double CDF(const double x) const;
		// compute iCDF(p) of N(mean, stdev)
		double iCDF(const double p) const;

		/*
		the following functions are associated with a normal distribution
		with mean =0 and stdev = 1
		*/

		static double pdf(const double x);

		static double moro_cdf(const double x);
		static double cody_cdf(const double x);
		static double intfrac_cody_cdf(const double x, double *intpart, double *fracpart);
		static double density_cody_cdf(const double x, double *density);

		static double moro_icdf(const double x);
		static double hp_moro_icdf(const double x);
		static double acklam_icdf(const double p);

		static constexpr double OneOverRootTwoPi = 0.398942280401433;
		static constexpr double OneOverRootTwo = 0.7071067811865;

	private:
		const double mean, stdev;
		double alpha;
		double beta;
		



	};
