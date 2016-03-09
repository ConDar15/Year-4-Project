#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#include "taylor_inv_trig.h"
#include "trig_utilities.h"
#include "utilities.h"

#define INIT_CONSTANTS in = fopen(PI_INFILE, "r"); \
					   mpfr_init(MPFR_PI); \
					   mpfr_inp_str(MPFR_PI, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   mpfr_init(MPFR_HALF_PI); \
					   mpfr_div_ui(MPFR_HALF_PI, MPFR_PI, 2, MPFR_RNDN);

mpfr_t MPFR_PI, MPFR_HALF_PI;

double taylor_asin(double x, unsigned int N)
{
	assert(x >= -1 && x <= 1);
	
	if(x == 1)
		return HALF_PI;
	else if(x == -1)
		return -1 * HALF_PI;
	
	double s = x, x_2 = x*x, t = x;

	for(int n = 1; n < N; n++)
	{
		t *= 2*n*(2*n - 1)*x_2;
		t /= 4*n*n;
		s += t/(2*n+1);
	}
	return s;
}

void mpfr_taylor_asin(mpfr_t R, mpfr_t x, unsigned int N)
{
	assert(mpfr_cmp_si(x, -1) >= 0 && mpfr_cmp_si(x, 1) <= 0);

	if(mpfr_cmp_si(x, 1) == 0)
		mpfr_set(R, MPFR_HALF_PI, MPFR_RNDN);
	else if(mpfr_cmp_si(x, -1) == 0)
	{
		mpfr_set(R, MPFR_HALF_PI, MPFR_RNDN);
		mpfr_neg(R, R, MPFR_RNDN);
	}
	else
	{
		mpfr_t x_2, t, T;
		mpfr_init(x_2);
		mpfr_mul(x_2, x, x, MPFR_RNDN);
		mpfr_init_set(t, x, MPFR_RNDN);
		mpfr_init(T);
		mpfr_set(R, x, MPFR_RNDN);
		for(int n = 1; n < N; n++)
		{
			mpfr_mul_ui(t, t, 2*n*(2*n - 1), MPFR_RNDN);
			mpfr_mul(t, t, x_2, MPFR_RNDN);
			mpfr_div_ui(t, t, 4*n*n, MPFR_RNDN);
			mpfr_div_ui(T, t, 2*n + 1, MPFR_RNDN);
			mpfr_add(R, R, T, MPFR_RNDN);
		}
	}
}

double taylor_acos(double x, unsigned int N)
{
	return HALF_PI - taylor_asin(x, N);
}

void mpfr_taylor_acos(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_taylor_asin(R, x, N);
	mpfr_sub(R, MPFR_HALF_PI, R, MPFR_RNDN);
}

double taylor_atan_bounded(double x, unsigned int N)
{
	assert(x >= 0 && x <= 1);
	double t = 0, x_2 = x*x, y = x;
	for(int n = 0; n < N; n++)
	{
		t += y/(2*(n++) + 1);
		y *= x_2;
		t -= y/(2*n + 1);
		y *= x_2;
	}
	return t;
}

void mpfr_taylor_atan_bounded(mpfr_t R, mpfr_t x, unsigned int N)
{
	assert(mpfr_cmp_ui(x, 0) >= 0 && mpfr_cmp_ui(x, 1) <= 1);
	mpfr_t t, x_2, y, a;
	mpfr_init_set(y, x, MPFR_RNDN);
	mpfr_init_set_ui(t, 0, MPFR_RNDN);
	mpfr_init(x_2);
	mpfr_mul(x_2, x, x, MPFR_RNDN);
	mpfr_init(a);
	for(int n = 0; n < N; n++)
	{
		mpfr_div_ui(a, y, 2*(n++) + 1, MPFR_RNDN);
		mpfr_add(t, t, a, MPFR_RNDN);
		mpfr_mul(y, y, x_2, MPFR_RNDN);
		mpfr_div_ui(a, y, 2*n + 1, MPFR_RNDN);
		mpfr_sub(t, t, a, MPFR_RNDN);
		mpfr_mul(y, y, x_2, MPFR_RNDN);
	}
}

double taylor_atan(double x, unsigned int N)
{
	if(x < 0)
		return taylor_atan(-x, N);
	
	if(x > 1)
		return HALF_PI/2 + taylor_atan_bounded((x - 1)/(x + 1), N);
	
	return taylor_atan_bounded(x, N);
}

void mpfr_taylor_atan(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_t y, pi_4, z;
	mpfr_init(y);
	mpfr_init(pi_4);
	mpfr_init(z);

	if(mpfr_cmp_ui(x, 0) < 0)
	{
		mpfr_neg(y, x, MPFR_RNDN);
		mpfr_taylor_atan(R, y, N);
	}
	else if(mpfr_cmp_ui(x, 1) > 0)
	{
		mpfr_div_ui(pi_4, MPFR_HALF_PI, 2, MPFR_RNDN);
		mpfr_add_ui(z, x, 1, MPFR_RNDN);
		mpfr_sub_ui(y, x, 1, MPFR_RNDN);
		mpfr_div(y, y, z, MPFR_RNDN);
		mpfr_taylor_atan_bounded(R, y, N);
		mpfr_add(R, R, pi_4, MPFR_RNDN);
	}
	else
	{
		mpfr_taylor_atan_bounded(R, y, N);
	}
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double x, y;
	unsigned int n, p, D;
	mpfr_t R, X;
	char sf[50];
	FILE *in;
	
	if(argc > 1)
	{
		switch(argv[1][0])
		{
			case 'a':
				if(argc == 5 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("arcCos(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_acos(x, n));
				else
					printf("Usage: %s a <x=value for arcCos(x)> <n> "
						   "<D=Number of digits to display>\n",
							argv[0]);
				break;
			 
			case 'b':
				if(argc == 5 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("arcSin(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_asin(x, n));
				else
					printf("Usage: %s b <x=value for arcSin(x)> <n> "
						   "<D=Number of digits to display>\n",
							argv[0]);
				break;

			case 'c':
				if(argc == 5 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Tan(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_atan(x, n));
				else
					printf("Usage: %s a <x=value for arcTan(x)> <n> "
						   "<D=Number of digits to display>\n",
							argv[0]);
				break;
			
			case 'd':
				if(argc == 6 && 
				   sscanf(argv[3], "%u", &D) == 1 &&
				   sscanf(argv[4], "%u", &n) == 1 &&
				   sscanf(argv[5], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS
					if (mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "arcCos(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_acos(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for arcCos(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for arcCos(x)> "
						   "<D=Number of digits to display> "
						   "<n> <p=bits of precision to use>\n",
						   argv[0]);
				break;
			
			case 'e':
				if(argc == 6 && 
				   sscanf(argv[3], "%u", &D) == 1 &&
				   sscanf(argv[4], "%u", &n) == 1 &&
				   sscanf(argv[5], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS
					if (mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "arcSin(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_asin(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for arcSin(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for arcSin(x)> "
						   "<D=Number of digits to display> "
						   "<n> <p=bits of precision to use>\n",
						   argv[0]);
				break;
			
			case 'f':
				if(argc == 6 && 
				   sscanf(argv[3], "%u", &D) == 1 &&
				   sscanf(argv[4], "%u", &n) == 1 &&
				   sscanf(argv[5], "%u", &p) == 1)
				{
					mpfr_set_default_prec(p);
					INIT_CONSTANTS
					if (mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN) == 0)
					{
						mpfr_init(R);

						sprintf(sf, "arcTan(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_atan(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for arcTan(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for arcTan(x)> "
						   "<D=Number of digits to display> "
						   "<n> <p=bits of precision to use>\n",
						   argv[0]);
				break;

			default:
				printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
}
#endif

