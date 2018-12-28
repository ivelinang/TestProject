#define BOOST_TEST_MODULE mytests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/math/distributions/normal.hpp>

#include "NormalDistribution1D.h"
#include <iostream>


BOOST_AUTO_TEST_CASE(myTestCase)
{
  BOOST_CHECK(1 == 1);
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(NormalDistributionTest)
{
	NormalDistribution1D myNorm(0, 1);
	double p = 4.5;
	double res = myNorm.CDF(p);
	boost::math::normal norm;
	double boost_res = boost::math::cdf(norm, p);
	BOOST_CHECK_CLOSE(boost_res, res, 0.0001);

	double inv_res = myNorm.iCDF(res);
	double boost_inv_res = boost::math::quantile(norm, res);

	BOOST_CHECK_CLOSE(inv_res, boost_inv_res, 0.0001);
}


BOOST_AUTO_TEST_CASE(InverseCDFTest)
{
	NormalDistribution1D myNorm(0, 1);
	boost::math::normal norm;
	double p = 0.65;
	double inv_res = myNorm.iCDF(p);
	double boost_inv_res = boost::math::quantile(norm, p);

	BOOST_CHECK_CLOSE(inv_res, boost_inv_res, 0.0001);
}

