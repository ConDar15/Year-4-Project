#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <assert.h>
#include <math.h>

#include "utilities.h"
#include "exact_root.h"
#include "newton_sqrt.h"

#define INIT_CONSTANTS mpfr_init_set_d(MPFR_HALF, 0.5, MPFR_RNDN); \
					   in = fopen(ROOT_2_INFILE, "r"); \
					   mpfr_init(MPFR_ROOT_2); \
					   mpfr_inp_str(MPFR_ROOT_2, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   in = fopen(ROOT_2_INV_INFILE, "r"); \
					   mpfr_init(MPFR_ROOT_2_INV); \
					   mpfr_inp_str(MPFR_ROOT_2_INV, in, 10, MPFR_RNDN); \
					   fclose(in);

mpfr_t MPFR_ROOT_2, MPFR_ROOT_2_INV, MPFR_HALF;

double newton_sqrt_v1(double N, double T)
{
	assert(N >= 0);
	assert(T >= 0);

	double x, px, d;

	//sets the initial guess for x
	x = N > 1 ? N : 1;
	//initial iterative error set to a high value to ensure it is larger
	//	than the given tolerance, provided T is a reasonable tolerance.
	d = 1000000;

	//Continues while the current iterative error is more than the
	//	provided tolerance.
	while(d > T)
	{
		//Updates store of the previously approximation
		px = x;
		//Calculates the next approximation
		x = 0.5 * (x + N/x);
		//Updates the iterative error
		d = fabs(x - px);
	}
	return x;
}

double newton_sqrt_v2(double N, double T)
{
	assert(N >= 0);
	assert(T >= 0);

	double x, px, d;

	//Sets the initial guess as explained in the report section on
	//	Newton Square Root approximations
	if(N >= 4)
		x = uint_sqrt((unsigned long) N);
	else
		x = N > 1 ? N : 1;

	//initial iterative error set to a high value to ensure it is larger
	//	than the given tolerance, provided T is a reasonable tolerance.
	d = 1000000;

	//Continues while the current iterative error is more than the
	//	provided tolerance.
	while(d > T)
	{
		//Updates store of the previously approximation
		px = x;
		//Calculates the next approximation
		x = 0.5 * (x + N/x);
		//Updates the iterative error
		d = fabs(x - px);
	}

	return x;
}

double newton_sqrt_v3(double N, double T)
{
	assert(N >= 0);
	assert(T >= 0);

	int e;
	double x, px, d;

	//frexp ensures 1/2 <= N < 1 and N*2^e = input N	
	N = frexp(N, &e);

	//Initial approximation is now a fixed constant
	x = 1;
	//initial iterative error set to a high value to ensure it is larger
	//	than the given tolerance, provided T is a reasonable tolerance.
	d = 1000000;

	//Continues while the current iterative error is more than the
	//	provided tolerance.
	while(d > T)
	{
		//Updates store of the previously approximation
		px = x;
		//Calculates the next approximation
		x = 0.5 * (x + N/x);
		//Updates the iterative error
		d = fabs(x - px);
	}

	//If e is odd then the result must be multiplied by either sqrt(2) or
	//	1/sqrt(2) to give the correct approximation
	if(e%2)
		x *= e > 0 ? ROOT_2 : ROOT_2_INV;
	//return x*2^(e/2)
	return ldexp(x, e / 2);
}

double newton_sqrt_v3_it(double N, unsigned int I)
{
	assert(N >= 0);

	int e;
	double x;
	
	N = frexp(N, &e);

	x = 1;

	for(int i = 0; i < I; ++i)
		x = 0.5 * (x + N/x);

	if(e%2)
		x *= e > 0 ? ROOT_2 : ROOT_2_INV;
	return ldexp(x, e / 2);
}


void mpfr_newton_sqrt_v3(mpfr_t R, mpfr_t N, mpfr_t T)
{
	mpfr_t x, px, d, t, n;
	mpfr_exp_t e;

	mpfr_init(n);
	mpfr_frexp(&e, n, N, MPFR_RNDN);

	mpfr_init_set_ui(x, 1, MPFR_RNDN);
	mpfr_init(px);
	mpfr_init_set_ui(d, 1000000, MPFR_RNDN);
	mpfr_init(t);

	while(mpfr_cmp(d, T) > 0)
	{
		mpfr_set(px, x, MPFR_RNDN);
		mpfr_div(t, n, x, MPFR_RNDN);
		mpfr_add(x, x, t, MPFR_RNDN);
		mpfr_mul(x, MPFR_HALF, x, MPFR_RNDN);
		mpfr_sub(d, x, px, MPFR_RNDN);
		mpfr_abs(d, d, MPFR_RNDN);
	}

	if(e%2)
		if(e > 0)
			mpfr_mul(x, MPFR_ROOT_2, x, MPFR_RNDN);
		else
			mpfr_mul(x, MPFR_ROOT_2_INV, x, MPFR_RNDN);
	mpfr_mul_2si(R, x, e/2, MPFR_RNDN);
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double N, T;
	unsigned int n, D, p;
	mpfr_t Nr, Tr, R;
	int c;
	char sf[50];
	FILE *in;

	if(argc==1)
	{
		printf("Usage: %s [a/b/c/d/e] <Arguments>\n", argv[0]);
		exit(1);
	}

	switch(argv[1][0])
	{
		case 'a':
			if (argc == 5 &&
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u" , &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, newton_sqrt_v1(N, T));
			else
				printf("Usage: %s a <N=Value to sqrt> "
					   "<T=Tolerance> <D=Number of digits to display>\n",
					   argv[0]);
			break;
		
		case 'b':
			if (argc == 5 &&
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u" , &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, newton_sqrt_v2(N, T));
			else
				printf("Usage: %s b <N=Value to sqrt> "
					   "<T=Tolerance> <D=Number of digits to display>\n",
					   argv[0]);
			break;

		case 'c':
			if (argc == 5 &&
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u" , &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, newton_sqrt_v3(N, T));
			else
				printf("Usage: %s c <N=Value to sqrt> "
					   "<T=Tolerance> <D=Number of digits to display>\n",
					   argv[0]);
			break;

		case 'd':
			if (argc == 5 &&
					sscanf(argv[3], "%u", &D) == 1 &&
					sscanf(argv[4], "%u" , &p) == 1)
			{
				mpfr_set_default_prec(p);
				INIT_CONSTANTS
							
				if (mpfr_init_set_str(Nr, argv[2], 10, MPFR_RNDN) == 0)
				{
					mpfr_init(R);

					mpfr_digits_to_tolerance(D, Tr);
					
					sprintf(sf, "sqrt(%%.%uRNf) =~\n\t%%.%uRNf\n", d(D), D);
					
					mpfr_newton_sqrt_v3(R, Nr, Tr);
					mpfr_printf(sf, Nr, R);
				}
				else
					printf("Usage: %s d <N=Value to sqrt> "
						   "<D=Number of digitsto calculate to> "
						   "<p=bits of precision>\n", argv[0]);
			}
			else
				printf("Usage: %s d <N=Value to sqrt> "
					   "<D=Number of digitsto calculate to> "
					   "<p=bits of precision>\n", argv[0]);
			break;
		
		case 'e':
			if (argc == 5 &&
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%u",  &p) == 1 &&
					sscanf(argv[4], "%u" , &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, 
												 newton_sqrt_v3_it(N, p));
			else
				printf("Usage: %s e <N=Value to sqrt> "
					   "<I=Iterations> <D=Number of digits to display>\n",
					   argv[0]);
			break;
		default:
			printf("Usage: %s [a/b/c/d/e] <Arguments>\n", argv[0]);
	}
}
#endif
