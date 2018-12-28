#pragma once

#include <vector>

double wave(double x, double y);

void wave_add(const double x, const double y, const double d_res, const double d_dx, const double d_dy);

double g(int i, double x1, double x2);

void h(const std::vector<double>& x, double& r1, double& r2);

double f(const std::vector<double>& x);


void h_aad(const std::vector<double>& x, double d_dr1, double d_dr2, std::vector<double>& d_dx);

void g_aad(int i, double x1, double x2, double d_dresult, double& d_dx1, double& d_dx2);

void f_aad(const std::vector<double>& x, double d_dresult, std::vector<double>& d_dx);


// f(x1, x2) = log(x1*x2)*x1


double f_m(double x1, double x2);

double g_m(double x1, double x2);

double f_m_ad(double x1, double x2, double dx1_d, double dx2_d);

double g_m_ad(double x1, double x2, double dx1_d, double dx2_d);

void f_m_aad(double x1, double x2, double d_dresult, double& d_dx1, double& d_dx2);

double f_m_aad2(double x1, double x2, double d_f, double d_dx1, double d_dx2);

void wave_aad_2(const double x, const double y, const double d_dresult, double& d_dx, double& d_dy);
