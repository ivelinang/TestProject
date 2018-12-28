#pragma once

#include <iostream>
#include <vector>

struct cashflow
{
	double notional_;
	//date date_;
	double discount_;
};

struct cashflows_v2
{
	std::vector<double> notional_;
	std::vector<double> discount_;
};

double pv_cashflows(const std::vector<cashflow>& flows);

double pv_cashflows_v2(const cashflows_v2 flows);
