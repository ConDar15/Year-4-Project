#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

#include "int_exp.h"
#include "log_exp_utilities.h"
#include "utilities.h"
#include "taylor_exp_log.h"

#define INIT_CONSTANTS in = fopen(NAT_LOG_2_INFILE, "r"); \
					   mpfr_init(MPFR_NAT_LOG_2); \
					   mpfr_inp_str(MPFR_NAT_LOG_2, in, 10, MPFR_RNDN); \
					   fclose(in);

mpfr_t MPFR_NAT_LOG_2;

double naive_exp(double x, unsigned int n)
{
	double z = 1 + x/n;
	return z > 0 ? squaring_int_exp(z, n)
				 : n%2 ? -squaring_int_exp(-z, n)
				 	   : squaring_int_exp(-z, n);
}

void mpfr_naive_exp(mpfr_t R, mpfr_t x, mpz_t n)
{
	assert(mpz_cmp_ui(n, 0) >= 0);
	
	mpfr_t z;
	
	mpfr_init_set(z, x, MPFR_RNDN);
	mpfr_div_z(z, z, n, MPFR_RNDN);
	mpfr_add_ui(z, z, 1, MPFR_RNDN);
	
	if(mpfr_cmp_ui(z, 0) > 0)
		mpfr_squaring_int_exp(R, z, n);
	else
	{
		mpfr_neg(z, z, MPFR_RNDN);
		mpfr_squaring_int_exp(R, z, n);
		
		if(mpz_odd_p(n))
			mpfr_neg(R, R, MPFR_RNDN);
	}
}

double taylor_exp(double x, unsigned int n)
{	
	double t, z;
	t = 1;
	z = 1;
	
	for(int k = 1; k < n; ++k)
	{
		t *= x;
		t /= k;
		z += t;
	}
	
	return z;
}

void mpfr_taylor_exp(mpfr_t R, mpfr_t x, mpz_t n)
{
	assert(mpz_cmp_ui(n, 0) >= 0);
	
	mpfr_t t;
	mpz_t k;
	
	mpfr_init_set_ui(t, 1, MPFR_RNDN);
	mpfr_set_ui(R, 1, MPFR_RNDN);
	
	for(mpz_init_set_ui(k, 1); mpz_cmp(k, n) < 0; mpz_add_ui(k, k, 1))
	{
		mpfr_mul(t, t, x, MPFR_RNDN);
		mpfr_div_z(t, t, k, MPFR_RNDN);
		mpfr_add(R, R, t, MPFR_RNDN);
	}
}

double taylor_nat_log(double x, unsigned int n)
{
	assert(n > 0);
	assert(x > 0);

	double a, t, z;
	int b;

	a = frexp(x, &b);
	a = 1 - a;
	
	t = a;
	z = a;
	
	for(int k = 2; k < n; ++k)
	{
		t *= a;
		z += t/k;
	}

	return b*NAT_LOG_2 - z;
}

void mpfr_taylor_nat_log(mpfr_t R, mpfr_t x, mpz_t n)
{
	assert(mpz_cmp_ui(n, 0) > 0);
	assert(mpfr_cmp_ui(x, 0) > 0);

	mpfr_t a, t, tt;
	mpfr_exp_t b;
	mpz_t k;
	unsigned int f = 1000, F = 1000;

	mpfr_init(a);
	mpfr_frexp(&b, a, x, MPFR_RNDN);
	mpfr_ui_sub(a, 1, a, MPFR_RNDN);
	
	mpfr_init_set(t, a, MPFR_RNDN);
	mpfr_init(tt);
	mpfr_set(R, a, MPFR_RNDN);

	for(mpz_init_set_ui(k, 2); mpz_cmp(k, n) < 0; mpz_add_ui(k, k, 1))
	{
		mpfr_mul(t, t, a, MPFR_RNDN);
		mpfr_div_z(tt, t, k, MPFR_RNDN);
		mpfr_add(R, R, tt, MPFR_RNDN);
	}

	mpfr_mul_si(a, MPFR_NAT_LOG_2, b, MPFR_RNDN);
	mpfr_sub(R, a, R, MPFR_RNDN);
}

double taylor_log(double x, double y, unsigned int n)
{
	assert(x > 0);
	assert(y > 0);
	assert(n > 0);

	return taylor_nat_log(y, n)/taylor_nat_log(x, n);
}

double taylor_pow(double x, double y, unsigned int n)
{
	assert(x > 0);
	assert(n > 0);

	return taylor_exp(y*taylor_nat_log(x, n), n);
}

void mpfr_taylor_log(mpfr_t R, mpfr_t x, mpfr_t y, mpz_t n)
{
	assert(mpfr_cmp_ui(x, 0) > 0);
	assert(mpfr_cmp_ui(y, 0) > 0);
	assert(mpz_cmp_ui(n, 0)  >= 0);

	mpfr_t A;
	mpfr_init(A);
	
	mpfr_taylor_nat_log(A, x, n);
	mpfr_taylor_nat_log(R, y, n);
	
	mpfr_div(R, R, A, MPFR_RNDN);
}

void mpfr_taylor_pow(mpfr_t R, mpfr_t x, mpfr_t y, mpz_t n)
{
	assert(mpfr_cmp_ui(x, 0) > 0);
	assert(mpz_cmp_ui(n, 0)  > 0);

	mpfr_t A;
	mpfr_init(A);

	mpfr_taylor_nat_log(A, x, n);
	mpfr_mul(A, y, A, MPFR_RNDN);
	mpfr_taylor_exp(R, A, n);
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
					printf("exp(%.*lf) ~= %.*lf (naive)\n", 
						d(D), x, D, naive_exp(x, n));
				else
					printf("Usage: %s a <x=Value for exp(x)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;

			case 'b':
				if(argc == 5 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%u",  &n) == 1 &&
					sscanf(argv[4], "%u",  &D) == 1)
					printf("exp(%.*lf) ~= %.*lf\n", 
						d(D), x, D, taylor_exp(x, n));
				else
					printf("Usage: %s b <x=Value for exp(x)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;
				
			case 'c':
				if(argc == 5 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%u",  &n) == 1 &&
					sscanf(argv[4], "%u",  &D) == 1)
					printf("ln(%.*lf) ~= %.*lf\n", 
						d(D), x, D, taylor_nat_log(x, n));
				else
					printf("Usage: %s c <x=Value for ln(x)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;

			case 'd':
				if(argc == 6 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%lf", &y) == 1 &&
					sscanf(argv[4], "%u",  &n) == 1 &&
					sscanf(argv[5], "%u",  &D) == 1)
					printf("pow(%.*lf, %.*lf) ~= %.*lf\n", 
						d(D), x, d(D), y, D, taylor_pow(x, y, n));
				else
					printf("Usage: %s d <x=Value for pow(x,y)> "
						   "<y=Value for pow(x,y)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;

			case 'e':
				if(argc == 6 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%lf", &y) == 1 &&
					sscanf(argv[4], "%u",  &n) == 1 &&
					sscanf(argv[5], "%u",  &D) == 1)
					printf("log(%.*lf, %.*lf) ~= %.*lf\n", 
						d(D), x, d(D), y, D, taylor_log(x, y, n));
				else
					printf("Usage: %s e <x=Value for log(x, y)> "
						   "<y=Value for log(x, y)> <n> "
						   "<D=digits to display>\n", argv[0]);
				break;

			case 'f':
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

						sprintf(sf, "exp(%%.%uRNf) =~\t(naive)"
									 "\n\t%%.%uRNf\n", d(D), D);
						
						mpfr_naive_exp(R, X, N);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s f <X=value for exp(X)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s f <X=value for exp(X)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
			
			case 'g':
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

						sprintf(sf, "exp(%%.%uRNf) =~"
									 "\n\t%%.%uRNf\n", d(D), D);
						
						mpfr_taylor_exp(R, X, N);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s g <X=value for exp(X)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s g <X=value for exp(X)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
				
			case 'h':
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
						
						sprintf(sf, "ln(%%.%uRNf) =~"
									 "\n\t%%.%uRNf\n", d(D), D);
						
						mpfr_taylor_nat_log(R, X, N);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s h <X=value for ln(X)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s h <X=value for ln(X)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
				
			case 'i':
				if(argc == 7 &&
					sscanf(argv[5], "%u", &D) == 1 &&
					sscanf(argv[6], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS

					if(mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0 &&
						mpfr_init_set_str(Y, argv[3],10, MPFR_RNDN) == 0 &&
						mpz_init_set_str(N, argv[4], 10) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "pow(%%.%uRNf, %%.%uRNf) =~"
									 "\n\t%%.%uRNf\n", d(D), d(D), D);
						
						mpfr_taylor_pow(R, X, Y, N);
						mpfr_printf(sf, X, Y, R);
					}
					else
						printf("Usage: %s i <X=value for pow(X,Y)> "
							   "<Y=value for pow(X,Y)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s i <X=value for pow(X,Y)> "
						   "<Y=value for pow(X,Y)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
				
			case 'j':
				if(argc == 7 &&
					sscanf(argv[5], "%u", &D) == 1 &&
					sscanf(argv[6], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS

					if(mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0 &&
						mpfr_init_set_str(Y, argv[3],10, MPFR_RNDN) == 0 &&
						mpz_init_set_str(N, argv[4], 10) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "log(%%.%uRNf, %%.%uRNf) =~"
									 "\n\t%%.%uRNf\n", d(D), d(D), D);
						
						mpfr_taylor_log(R, X, Y, N);
						mpfr_printf(sf, X, Y, R);
					}
					else
						printf("Usage: %s j <X=value for log(X,Y)> "
							   "<Y=value for log(X,Y)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s j <X=value for log(X,Y)> "
						   "<Y=value for log(X,Y)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
			
			default:
				printf("Usage: %s <a/b/c/d/e/f/g/h/i/j> <arguments>\n", 
					argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c/d/e/f/g/h/i/j> <arguments>\n", argv[0]);		
}
#endif
