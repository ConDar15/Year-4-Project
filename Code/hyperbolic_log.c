#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

#include "log_exp_utilities.h"
#include "utilities.h"
#include "hyperbolic_log.h"

#define INIT_CONSTANTS in = fopen(NAT_LOG_2_INFILE, "r"); \
					   mpfr_init(MPFR_NAT_LOG_2); \
					   mpfr_inp_str(MPFR_NAT_LOG_2, in, 10, MPFR_RNDN); \
					   fclose(in);

mpfr_t MPFR_NAT_LOG_2;

double hyperbolic_nat_log(double x, unsigned int n)
{
	//Ensuers that x is positive
	assert(x > 0);

	double t, y, z, a;
	unsigned int d;
	int b;

	//Ensures 1/2 <= a < 1 and a*2^b == x
	a = frexp(x, &b);
	
	//Sets the initial values to be used
	t = (a-1)/(a+1);
	d = 1;
	y = t*t;
	z = t;
	
	//Main loop
	for(unsigned int k = 1; k <= n; ++k)
	{
		//Performs a single Taylor Series update step
		t *= y;
		d += 2;
		z += t/d;
	}
	
	return b*NAT_LOG_2 + 2*z;
}

void mpfr_hyperbolic_nat_log(mpfr_t R, mpfr_t x, mpz_t n)
{
	assert(mpfr_cmp_ui(x, 0) > 0);
	assert(mpz_cmp_ui(n, 0) >= 0);

	mpfr_t t, y, z, tmp, a;
	mpfr_exp_t b;
	mpz_t d, k;
	
	mpfr_init(a);
	mpfr_init(t);
	mpfr_init(tmp);
	mpfr_init(y);

	mpfr_frexp(&b, a, x, MPFR_RNDN);
	
	mpfr_sub_ui(t, a, 1, MPFR_RNDN);
	mpfr_add_ui(tmp, a, 1, MPFR_RNDN);
	mpfr_div(t, t, tmp, MPFR_RNDN);
	mpfr_mul(y, t, t, MPFR_RNDN);
	
	mpfr_init_set(z, t, MPFR_RNDN);
	mpz_init_set_ui(d, 1);
	
	for(mpz_init_set_ui(k, 1); mpz_cmp(k, n) <= 0; mpz_add_ui(k, k, 1))
	{
		mpfr_mul(t, t, y, MPFR_RNDN);
		mpz_add_ui(d, d, 2);
		mpfr_div_z(tmp, t, d, MPFR_RNDN);
		mpfr_add(z, z, tmp, MPFR_RNDN);
	}
	
	mpfr_mul_ui(R, z, 2, MPFR_RNDN);
	mpfr_mul_si(tmp, MPFR_NAT_LOG_2, b, MPFR_RNDN);
	mpfr_add(R, R, tmp, MPFR_RNDN);
}

double hyperbolic_log(double x, double y, unsigned int n)
{
	assert(x > 0);
	assert(y > 0);

	return hyperbolic_nat_log(y, n)/hyperbolic_nat_log(x, n);
}

void mpfr_hyperbolic_log(mpfr_t R, mpfr_t x, mpfr_t y, mpz_t n)
{
	assert(mpfr_cmp_ui(x, 0) > 0);
	assert(mpfr_cmp_ui(y, 0) > 0);
	assert(mpz_cmp_ui(n, 0) >= 0);

	mpfr_t A;
	mpfr_init(A);

	mpfr_hyperbolic_nat_log(R, y, n);
	mpfr_hyperbolic_nat_log(A, x, n);
	
	mpfr_div(R, R, A, MPFR_RNDN);
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double x, y;
	unsigned int n, D, p;
	mpfr_t X, Y, R;
	mpz_t N;
	char sf[50];
	FILE *in;

	if(argc > 1)
	{
		switch(argv[1][0])
		{
			case 'a':
				if(argc == 5 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%u",  &n) == 1 &&
					sscanf(argv[4], "%u",  &D) == 1)
					printf("ln(%.*lf) ~= %.*lf\n",
						d(D), x, D, hyperbolic_nat_log(x, n));
				else
					printf("Usage: %s a <x=Value for ln(x)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;
			
			case 'b':
				if(argc == 6 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%lf", &y) == 1 &&
					sscanf(argv[4], "%u",  &n) == 1 &&
					sscanf(argv[5], "%u",  &D) == 1)
					printf("log(%.*lf, %.*lf) ~= %.*lf\n",
						d(D), x, d(D), y, D, hyperbolic_log(x, y, n));
				else
					printf("Usage: %s b <x=Value for log(x,y)> "
						   "<y=Value for log(x,y)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;

			case 'c':
				if(argc == 6 &&
					sscanf(argv[4], "%u", &D) == 1 &&
					sscanf(argv[5], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS

					if(mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0 &&
						mpz_init_set_str(N, argv[3], 10) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "ln(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);

						mpfr_hyperbolic_nat_log(R, X, N);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s c <X=value for ln(X)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s c <X=value for ln(X)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;

			case 'd':
				if(argc == 7 &&
					sscanf(argv[5], "%u", &D) == 1 &&
					sscanf(argv[6], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS

					if(mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0 &&
					   mpfr_init_set_str(Y, argv[3], 10, MPFR_RNDN) == 0 &&
					   	mpz_init_set_str(N, argv[4], 10) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "ln(%%.%uRNf, %%.%uRNf) =~\n\t"
									"%%.%uRNf\n", d(D), d(D), D);

						mpfr_hyperbolic_log(R, X, Y, N);
						mpfr_printf(sf, X, Y, R);
					}
					else
						printf("Usage: %s d <X=Value for ln(X,Y)> "
							   "<Y=Value for ln(X,Y)> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s d <X=Value for ln(X, Y)> "
						   "<Y=Value for ln(X,Y)> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;

			default:
				printf("Usage: %s <a/b/c/d> <arguments>\n", argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c/d> <arguments>\n", argv[0]);
}
#endif
