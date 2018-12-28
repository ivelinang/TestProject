
#include <iostream>
#include <algorithm>
#include <iterator>


#include "cashflow.h"
#include "fibonacci.h"
#include "smart_pointer.h"
#include "check_avx.h"
#include "test_shared_pointer.h"
#include "automatic_differentation.h"

#include <boost/lambda/lambda.hpp>




void f()
{
	std::cout << "Hello";
	std::cin.get();
}

void test_cashflow()
{
	cashflow cf1;
	cashflow cf2;
	cf1.notional_ = 100;
	cf1.discount_ = 0.99;

	std::vector<cashflow> v_cf{ cashflow(), cashflow(), cashflow() };

	for (int i = 0; i < v_cf.size(); ++i)
	{
		v_cf[i].notional_ = 100 + i * 10;
		v_cf[i].discount_ = 0.90 + i/100.0;

	}

	double res = pv_cashflows(v_cf);

	std::cout << "res: "<<res;
	std::cin.get();

	cashflows_v2 cf_v2;
	cf_v2.notional_ = { 100, 110, 120 };
	cf_v2.discount_ = { 0.90, 0.91, 0.92 };

	double res2 = pv_cashflows_v2(cf_v2);

	std::cout << "res2: " << res2;
	std::cin.get();
}

void test_fib()
{
	int x = fib(5);
	int y = fib_ext(5);

	std::cout << "slow fib: " << x << "\n";
	std::cout << "fast fib: " << y;
	std::cin.get();
}

void test_smart_pointer()
{
	SmartPointer<std::vector<int>> sp(new std::vector<int>{ 1,2,3 });
	std::cout<<sp->size();
	std::cin.get();
	{
		SmartPointer<std::vector<int>> q = sp;
		std::cout<<q->size();
		std::cin.get();
	}
	std::cout<<sp->size();
	std::cin.get();

}

void test_ad()
{
	double res = f_m_ad(2, 1, 1, 0);
	std::cout << "df_dx1 at x1=1: " << res;
	std::cin.get();

	double res2 = f_m_aad2(2, 1, 1, 1, 0);
	std::cout << "df_dx1 at x1=1 (v2): " << res2;
	std::cin.get();

	const double PI = std::atan(1.0) * 4;
	double x = PI / 4.0;
	double y = PI / 4.0;

	double d_dx = 0;
	double d_dy = 0;

	wave_aad_2(x, y, 1, d_dx, d_dy);

	std::cout << "d_dx : " << d_dx;
	std::cout << "d_dy : " << d_dy;
	std::cin.get();


}

void boost_test() {
	using namespace boost::lambda;
	typedef std::istream_iterator<int> in;

	std::for_each(
		in(std::cin), in(), std::cout << (_1 * 3) << " ");
}

/*
int main() {
	boost_test();
	f();
	checkAVX();
	test_cashflow();
	test_fib();
	test_smart_pointer();
	test_shared_pointer();
	test_ad();

	std::cout << "before main exit";
	std::cin.get();
	return 0;
}
*/