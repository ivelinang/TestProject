#include "fibonacci.h"

int fib(int x)
{
	if (x == 0)
		return 0;
	if (x == 1)
		return 1;

	return fib(x - 1) + fib(x - 2);
}

int fib_ext(int x)
{
	return fib_int(1, 1, x);
}

int fib_int(int prev, int current, int cnter)
{
	if (cnter == 1)
		return current;
	int x = prev;
	int y = current;
	prev = y;
	current = x + y;
	return fib_int(prev, current, --cnter);
}