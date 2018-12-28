
#include "cashflow.h"

double pv_cashflows(const std::vector<cashflow>& flows)
{
	double res = 0.0;
	for (size_t i = 0, n = flows.size(); i < n; ++i)
		res += flows[i].notional_ * flows[i].discount_;

	return res;
}

double pv_cashflows_v2(const cashflows_v2 flows)
{
	double res = 0.0;
	for (size_t i = 0; i < flows.notional_.size(); ++i)
		res += flows.notional_[i] * flows.discount_[i];

	return res;
}