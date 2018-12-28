
#include "automatic_differentation.h"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <numeric>

double wave(double x, double y)
{
	const double w1 = sin(x);
	const double w2 = cos(y);
	return w1 * w2;
}

void wave_add(const double x, const double y, const double d_dres, double d_dx, double d_dy)
{
	const double w1 = sin(x);
	const double dw1_dx = cos(x);

	const double w2 = cos(y);
	const double dw2_dy = -sin(y);

	const double dresult_dw1 = w2;
	const double dresult_dw2 = w1;

	//backward sweep of the graph
	const auto d_dw1 = d_dres * dresult_dw1;
	const auto d_dw2 = d_dres * dresult_dw2;

	d_dx += d_dw1 * dw1_dx;
	d_dy += d_dw2 * dw2_dy;


}


double f(const std::vector<double>& x)
{
	assert(x.size() >= 2);
	double r1 = 0;
	double r2 = 0;

	h(x, r1, r2);
	const double gr1 = g(x.size(), r1, r2);
	const double gr2 = g(2, x[0], x[1]);
	const double sum = std::accumulate(x.begin(), x.end(), 0.0);

	return sum * gr1*gr2;
}

double g(int i, double x1, double x2)
{
	return 3;
}

void h(const std::vector<double>& x, double& r1, double& r2)
{
	int dfd = 2;
}

void h_aad(const std::vector<double>& x, double d_dr1, double d_dr2, std::vector<double>& d_dx)
{
	int dfdf = 2;
}

void g_aad(int i, double x1, double x2, double d_dresult, double& d_dx1, double& d_dx2)
{
	int dfdfd = 2;
}


void f_aad(const std::vector<double>& x, double d_dresult, std::vector<double>& d_dx)
{
	double r1 = 0; 
	double r2 = 0;
	h(x, r1, r2);
	const double gr1 = g(x.size(), r1, r2); 
	const double gr2 = g(2, x[0], x[1]);
	const double sum = std::accumulate(x.begin(), x.end(), 0.0);
	//result = sum*gr1*gr2;

	//partial derivatives
	const double dresult_dsum = gr1 * gr2;
	const double dresult_dgr1 = sum * gr2;
	const double dresult_dgr2 = sum * gr1;

	//adjoint accumulation
	double d_dsum = d_dresult * dresult_dsum;
	double d_dgr1 = d_dresult * dresult_dgr1;
	double d_dgr2 = d_dresult * dresult_dgr2;

	for (double& d_delement : d_dx)
		d_delement += d_dsum;

	g_aad(2, x[0], x[1], d_dgr2, d_dx[0], d_dx[1]);
	double d_dr1 = 0;
	double d_dr2 = 0;
	g_aad(x.size(), r1, r2, d_dgr1, d_dr1, d_dr2);
	h_aad(x, d_dr1, d_dr2, d_dx);
}


double f_m_ad(double x1, double x2, double dx1_d, double dx2_d)
{
	//const double w = log(x1*x2);
	const double w = g_m(x1, x2);

	//const double dw_dx1 = 1 / x1;
	const double dw_dx1 = g_m_ad(x1, x2, 1, 0);
	//const double dw_dx2 = 1 / x2;
	const double dw_dx2 = g_m_ad(x1, x2, 0, 1);

	const double dresult_dw = x1;
	const double dresult_dx1 = w;

	return dresult_dw * (dw_dx1 * dx1_d + dw_dx2 * dx2_d) + dresult_dx1 * dx1_d;
}

double f_m(double x1, double x2)
{
	const double w = log(x1*x2);
	return x1 * w;
}

double g_m(double x1, double x2)
{
	return log(x1*x2);
}

double g_m_ad(double x1, double x2, double dx1_d, double dx2_d)
{
	const double dg_dx1 = 1 / x1;
	const double dg_dx2 = 1 / x2;

	return dg_dx1 * dx1_d + dg_dx2 * dx2_d;

}

void f_m_aad(double x1, double x2, double d_dresult, double& d_dx1, double& d_dx2)
{
	const double w = log(x1*x2);

	const double dw_dx1 = 1 / x1;
	const double dw_dx2 = 1 / x2;

	const double dresult_dw = x1;
	const double dresult_dx1 = w;

	d_dx1 += d_dresult * (dresult_dw * dw_dx1 + dresult_dx1);
	d_dx2 += d_dresult * dresult_dw * dw_dx2;

}


double f_m_aad2(double x1, double x2, double d_f, double d_dx1, double d_dx2)
{
	const double w = log(x1*x2);
	
	double df_dw = x1;
	double df_dx1 = w;

	double dw_dx1 = 1 / x1;
	double dw_dx2 = 1 / x2;

	

	double d_dw = d_f * df_dw;
	double d_x1 = d_f * df_dx1 + d_dw * dw_dx1;
	double d_x2 = d_dw * dw_dx2;

	return d_x1 * d_dx1 + d_dx2 * d_dx2;
}

void wave_aad_2(const double x, const double y, const double d_dresult, double& d_dx, double& d_dy)
{
	//wave(x,y) = sin(x)*cos(y)

	//forward step
	const double w1 = sin(x);
	const double w2 = cos(y);

	const double dw1_dx = cos(x);
	const double dw2_dy = -sin(y);

	const double dres_dw1 = w2;
	const double dres_dw2 = w1;

	//backward step
	const double d_dw1 = d_dresult * dres_dw1;
	const double d_dw2 = d_dresult * dres_dw2;

	d_dx += d_dw1 * dw1_dx;
	d_dy += d_dw2 * dw2_dy;
}