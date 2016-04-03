#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <assert.h>
#include <math.h>

#include "bisect_root.h"
#include "utilities.h"

#define INIT_CONSTANTS mpfr_init_set_ui(MPFR_ONE, 1, MPFR_RNDN); \
				       mpfr_init_set_d(MPFR_HALF, 0.5, MPFR_RNDN);

mpfr_t MPFR_ONE, MPFR_HALF;

double bisect_sqrt(double N, double T)
{
	assert(N >= 0);
	assert(T >= 0);

	int e;
	double a, b, x, f;

	//frexp finds a,b such that a*2^b = N and 1/2 <= a < 1
	N = frexp(N, &e);
	//Corrects for the case of odd exponents
	//	e%2 is true when e is odd
	if(e%2)
	{
		N /= 2;
		e += 1;
	}
	
	//Sets the initial values of a and b
	a = 0;
	b = 1;

	x = 0.5*(a + b);
	f = x*x - N;

	//fabs(f) > T is our approximation of 
	//  f != 0, by using the given tolerance
	while(fabs(f) > T && b - a > T)
	{
		//Update of the bounds a and b
		if (f < 0)
			a = x;
		else
			b = x;
		
		//Update of x and f
		x = 0.5*(a + b);
		f = x*x - N;
	}
	
	return ldexp(x, e / 2);
}

double bisect_sqrt_it(double N, unsigned int I)
{
	assert(N >= 0);

	int e;
	double a, b, x, f;

	//frexp finds a,b such that a*2^b = N and 1/2 <= a < 1
	N = frexp(N, &e);
	//Corrects for the case of odd exponents
	//	e%2 is true when e is odd
	if(e%2)
	{
		N /= 2;
		e += 1;
	}
	
	//Sets the initial values of a and b
	a = 0;
	b = 1;

	x = 0.5*(a + b);
	f = x*x - N;

	//fabs(f) > T is our approximation of 
	//  f != 0, by using the given tolerance
	for(int i = 0; i < I; ++i)
	{
		//Update of the bounds a and b
		if (f < 0)
			a = x;
		else
			b = x;
		
		//Update of x and f
		x = 0.5*(a + b);
		f = x*x - N;
	}
	
	return ldexp(x, e / 2);
}

double iPow(double x, unsigned int n)
{
	double r = 1;
	while(n--)
		r *= x;
	return r;
}

double bisect_nRoot(double N, double T, unsigned int n)
{
	assert(N >= 0);
	assert(T >= 0);
	//Ensures that none of the trivial cases are requested
	assert(n >= 2);

	//Runs the more optimal bisect_sqrt if n == 2
	if(n == 2)
		return bisect_sqrt(N, T);

	double a, b, x, f;
	
	//Sets the initial values of a and b
	a = 0;
	//This statement is equivalent to
	//  if(N<1) b=1; else b=N;
	b = N < 1 ? 1 : N;

	x = 0.5*(a + b);
	f = iPow(x, n) - N;
	
	//fabs(f) > T is our approximation of 
	//  f != 0, by using the given tolerance
	while(fabs(f) > T && b - a > T)
	{
		//Update of the bounds a and b
		if (f < 0)
			a = x;
		else
			b = x;
		
		//Update of x and f
		x = 0.5*(a + b);
		f = iPow(x, n) - N;
	}
	
	return x;
}

void mpfr_bisect_sqrt(mpfr_t R, mpfr_t N, mpfr_t T)
{
	if(mpfr_cmp_ui(N, 0) < 0)
	{
		fprintf(stderr, "The value to square root must be non-negative\n");
		exit(-1);
	}
	if(mpfr_cmp_ui(T, 0) < 0)
	{
		fprintf(stderr, "The tolerance must be non-negative\n");
		exit(-1);
	}
	
	mpfr_exp_t e;
	mpfr_t a, b, x, f, d, fab, n;

	mpfr_init(n);
	mpfr_frexp(&e, n, N, MPFR_RNDN);
	if(e%2)
	{
		mpfr_div_ui(n, n, 2, MPFR_RNDN);
		e += 1;
	}

	//Set a == 0
	mpfr_init_set_ui(a, 0, MPFR_RNDN);
	
	//Set b == 1
	mpfr_init_set_ui(b, 1, MPFR_RNDN);

	//Set x = (a + b)/2
	mpfr_init(x);
	mpfr_add(x, a, b, MPFR_RNDN);
	mpfr_mul(x, x, MPFR_HALF, MPFR_RNDN);
	
	//Set f = x^2 - N and fab = |f|
	mpfr_init(f);
	mpfr_init(fab);
	mpfr_mul(f, x, x, MPFR_RNDN);
	mpfr_sub(f, f, N, MPFR_RNDN);
	mpfr_abs(fab, f, MPFR_RNDN);

	//Set d = b - a
	mpfr_init(d);
	mpfr_sub(d, b, a, MPFR_RNDN);

	while(mpfr_cmp(fab, T) > 0 && mpfr_cmp(d, T) > 0)
	{
		//Update the bounds, a and b
		if(mpfr_cmp_ui(f, 0) < 0)
			mpfr_set(a, x, MPFR_RNDN);
		else
			mpfr_set(b, x, MPFR_RNDN);

		//Update x
		mpfr_add(x, a, b, MPFR_RNDN);
		mpfr_mul(x, x, MPFR_HALF, MPFR_RNDN);
		
		//Update f and fab
		mpfr_mul(f, x, x, MPFR_RNDN);
		mpfr_sub(f, f, n, MPFR_RNDN);
		mpfr_abs(fab, f, MPFR_RNDN);
	}

	printf("beep");
	mpfr_mul_2si(R, x, e/2, MPFR_RNDN);
}

void mpfr_bisect_nRoot(mpfr_t R, mpfr_t N, mpfr_t T, unsigned int n)
{
	if(mpfr_cmp_ui(N, 0) < 0)
	{
		fprintf(stderr, "The value to square root must be non-negative\n");
		exit(-1);
	}
	if(mpfr_cmp_ui(T, 0) < 0)
	{
		fprintf(stderr, "The tolerance must be non-negative\n");
		exit(-1);
	}
	assert(n >= 2);

	mpfr_t a, b, x, f, d, fab;

	//Set a == 0
	mpfr_init_set_ui(a, 0, MPFR_RNDN);
	
	//Set b = max{1, N}
	mpfr_init(b);
	mpfr_max(b, MPFR_ONE, N, MPFR_RNDN);

	//Set x = (a + b)/2
	mpfr_init(x);
	mpfr_add(x, a, b, MPFR_RNDN);
	mpfr_mul(x, x, MPFR_HALF, MPFR_RNDN);
	
	//Set f = x^2 - N
	mpfr_init(f);
	mpfr_init(fab);
	mpfr_pow_ui(f, x, n, MPFR_RNDN);
	mpfr_sub(f, f, N, MPFR_RNDN);
	mpfr_abs(fab, f, MPFR_RNDN);

	//Set d = b - a
	mpfr_init(d);
	mpfr_sub(d, b, a, MPFR_RNDN);

	while(mpfr_cmp(fab, T) > 0 && mpfr_cmp(d, T) > 0)
	{
		//Update the bounds, a and b
		if(mpfr_cmp_ui(f, 0) < 0)
			mpfr_set(a, x, MPFR_RNDN);
		else
			mpfr_set(b, x, MPFR_RNDN);

		//Update x
		mpfr_add(x, a, b, MPFR_RNDN);
		mpfr_mul(x, x, MPFR_HALF, MPFR_RNDN);
		
		//Update f
		mpfr_pow_ui(f, x, n, MPFR_RNDN);
		mpfr_sub(f, f, N, MPFR_RNDN);
		mpfr_abs(fab, f, MPFR_RNDN);
	}

	mpfr_set(R, x, MPFR_RNDN);
}

#ifdef COMPILE_MAIN
int main(int argc, char** argv)
{
	double N, T;
	unsigned int n, D, p;
	mpfr_t Nr, Tr, R;
	int c;
	char sf[50];

	if (argc == 1)
	{
		printf("Usage: %s [a/b/c/d/e] [arguments]\n", argv[0]);
		exit(1);
	}

	switch(argv[1][0])
	{
		case 'a':
			if (argc == 5 && 
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u", &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, bisect_sqrt(N, T));
			else
				printf("Usage: %s a <N=Value to sqrt> "
					   "<T=Tolerance> <D=Number of digits to display>\n", 
					   argv[0]);
			break;
		
		case 'b':
			if (argc == 6 && 
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u", &n) == 1 &&
					sscanf(argv[5], "%u", &D) == 1)
				printf("%u_Root(%.*lf) =~ %.*lf\n", 
						n, d(D), N, D, bisect_nRoot(N, T, n));
			else
				printf("Usage: %s b <N=Value to root> <T=Tolerance>"
				 	   "<n=nth Root> <D=Number of digits to display>\n", 
					   argv[0]);
			break;

		case 'c':
			if (argc == 5 &&
					sscanf(argv[3], "%u", &D) == 1 &&
					sscanf(argv[4], "%u", &p) == 1)
			{
				mpfr_set_default_prec(p);
				INIT_CONSTANTS

				if (mpfr_init_set_str(Nr, argv[2], 10, MPFR_RNDN) == 0)
				{
					mpfr_init(R);
					//Sets the tolerance to Tr = 10^-D
					mpfr_digits_to_tolerance(D, Tr);
					
					//Generates the required format string
					sprintf(sf, "sqrt(%%.%uRNf) =~\n\t%%.%uRNf\n", d(D), D);

					mpfr_bisect_sqrt(R, Nr, Tr);
					mpfr_printf(sf, Nr, R);
				}
				else
					printf("Usage: %s c <N=Value to sqrt> "
					       "<D=Number of digits to calculate to> "
						   "<p=bits of precision>\n", argv[0]);
			}
			else
				printf("Usage: %s c <N=Value to sqrt> "
				       "<D=Number of digits to calculate to> "
					   "<p=bits of precision>\n", argv[0]);
			break;

		case 'd':
			if (argc == 6 &&
					sscanf(argv[3], "%u", &D) == 1 &&
					sscanf(argv[4], "%u", &n) == 1 &&
					sscanf(argv[5], "%u", &p) == 1)
			{
				mpfr_set_default_prec(p);
				INIT_CONSTANTS

				if (mpfr_init_set_str(Nr, argv[2], 10, MPFR_RNDN) == 0)
				{
					mpfr_init(R);

					//Sets the tolerance to Tr = 10^-D
					mpfr_digits_to_tolerance(D, Tr);
					
					//Generates the required format string
					sprintf(sf, "%%u_Root(%%.%uRNf) =~\n\t%%.%uRNf\n", d(D), D);
					
					mpfr_bisect_nRoot(R, Nr, Tr, n);
					mpfr_printf(sf, n, Nr, R);
				}
				else
					printf("Usage: %s d <N=Value to root> "
					       "<D=Number of digits to calculate to> "
						   "<n=nth root> <p=bits of precision>\n", argv[0]);
			}
			else
				printf("Usage: %s d <N=Value to root> "
				       "<D=Number of digits to calculate to> " 
					   "<n=nth root> <p=bits of precision>\n", argv[0]);
			break;

		case 'e':
			if (argc == 5 && 
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%u", &p) == 1 &&
					sscanf(argv[4], "%u", &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, 
												 bisect_sqrt_it(N, p));
			else
				printf("Usage: %s a <N=Value to sqrt> "
					   "<I=iterartions> <D=Number of digits to display>\n", 
					   argv[0]);
			break;
		default:
			printf("Usage: %s [a/b/c/d/e] [arguments]", argv[0]);
	}
}
#endif
