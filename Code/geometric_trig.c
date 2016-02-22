#include <stdio.h>
#include <assert.h>
#include "trig_utilities.h"

double four_pow(unsigned int n)
{
	double x = 1;
	while(n--)
		x*=4;
	return x;
}

double geometric_cos(double x, unsigned int n)
{
	assert(PI_2 > x && x > 0);
	double S = (x*x)/four_pow(n);
	for(int i = 0; i < n; i++)
		S = S*(4-S);
	return 1 - S/2;
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double X, Y;
	unsigned int n;
	sscanf(argv[1], "%lf", &X);	
	sscanf(argv[2], "%u", &n);
	printf("Cos(%.*lf) ~= %.*lf\n", n+5, X, n+5, geometric_cos(X, n));
}
#endif
