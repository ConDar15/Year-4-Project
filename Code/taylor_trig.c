#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpfr.h>

#include "taylor_trig.h"
#include "trig_utilities.h"
#include "utilities.h"

#define INIT_CONSTANTS in = fopen(PI_INFILE, "r"); \
					   mpfr_init(MPFR_PI); \
					   mpfr_inp_str(MPFR_PI, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   mpfr_init(MPFR_TWO_PI); \
					   mpfr_init(MPFR_HALF_PI); \
					   mpfr_div_ui(MPFR_HALF_PI, MPFR_PI, 2, MPFR_RNDN); \
					   mpfr_mul_ui(MPFR_TWO_PI, MPFR_PI, 2, MPFR_RNDN);

mpfr_t MPFR_PI, MPFR_HALF_PI, MPFR_TWO_PI;

double taylor_cos_bounded(double x, unsigned int N)
{
	assert(x >= 0 && x <= HALF_PI);
	double c = 1, x_2 = x*x, a = 1, b = 1;
	for(int n = 1; n < N; n++)
	{
		a /= (2*n - 1)*(2*(n++));
		b *= x_2;
		c -= a*b;
		a /= (2*n - 1)*(2*(n));
		b *= x_2;
		c += a*b;
	}
	return c;
}

void mpfr_taylor_cos_bounded(mpfr_t R, mpfr_t x, unsigned int N)
{
	assert(mpfr_cmp_ui(x, 0) >= 0 && mpfr_cmp(x, MPFR_HALF_PI) <= 0);
	mpfr_t t, x_2;
	
	mpfr_init_set_ui(R, 1, MPFR_RNDN);
	mpfr_init_set_ui(t, 1, MPFR_RNDN);
	mpfr_init(x_2);
	mpfr_mul(x_2, x, x, MPFR_RNDN);

	for(int n = 1; n < N; n++)
	{
		mpfr_div_ui(t, t, (2*n-1)*(2*(n++)), MPFR_RNDN);
		mpfr_mul(t, t, x_2, MPFR_RNDN);
		mpfr_sub(R, R, t, MPFR_RNDN);
		mpfr_div_ui(t, t, (2*n-1)*(2*n), MPFR_RNDN);
		mpfr_mul(t, t, x_2, MPFR_RNDN);
		mpfr_add(R, R, t, MPFR_RNDN);
	}
}

double taylor_sin_bounded(double x, unsigned int N)
{
	assert(x >= 0 && x <= HALF_PI);
	double s = x, x_2 = x*x, a = 1, b = x;
	for(int n = 1; n < N; n++)
	{
		a /= (2*n + 1)*(2*(n++));
		b *= x_2;
		s -= a*b;
		a /= (2*n + 1)*(2*n);
		b *= x_2;
		s += a*b;
	}
	return s;
}

void mpfr_taylor_sin_bounded(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_printf("%.20RNF\n", x);
	assert(mpfr_cmp_ui(x, 0) >= 0 && mpfr_cmp(x, MPFR_HALF_PI) <= 0);
	mpfr_t t, x_2;
	
	mpfr_init_set(R, x, MPFR_RNDN);
	mpfr_init_set(t, x, MPFR_RNDN);
	mpfr_init(x_2);
	mpfr_mul(x_2, x, x, MPFR_RNDN);

	for(int n = 1; n < N; n++)
	{
		mpfr_div_ui(t, t, (2*n+1)*(2*(n++)), MPFR_RNDN);
		mpfr_mul(t, t, x_2, MPFR_RNDN);
		mpfr_sub(R, R, t, MPFR_RNDN);
		mpfr_div_ui(t, t, (2*n+1)*(2*n), MPFR_RNDN);
		mpfr_mul(t, t, x_2, MPFR_RNDN);
		mpfr_add(R, R, t, MPFR_RNDN);
	}
}

double taylor_cos(double x, unsigned int N)
{
	if(x >= 0)
	{
		while(x >= TWO_PI)
			x -= TWO_PI;

		if(x >= PI)
			if(x - PI >= HALF_PI)
				return taylor_cos_bounded(TWO_PI - x, N);
			else
				return -1  * taylor_cos_bounded(x - PI, N);
		else
			if(x >= HALF_PI)
				return -1 * taylor_cos_bounded(PI - x, N);
			else
				return taylor_cos_bounded(x, N);
	}
	
	return taylor_cos(-x, N);
}

void mpfr_taylor_cos(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_t y, t;
	mpfr_init_set(y, x, MPFR_RNDN);
	mpfr_init(t);
	if(mpfr_cmp_ui(y, 0) >= 0)
	{
		while(mpfr_cmp(y, MPFR_TWO_PI) >= 0)
			mpfr_sub(y, y, MPFR_TWO_PI, MPFR_RNDN);
		
		if(mpfr_cmp(y, MPFR_PI) >= 0)
		{
			mpfr_sub(t, y, MPFR_PI, MPFR_RNDN);
			if(mpfr_cmp(t, MPFR_HALF_PI) >= 0)
			{
				mpfr_sub(y, MPFR_TWO_PI, y, MPFR_RNDN);
				mpfr_taylor_cos_bounded(R, y, N);
			}
			else
			{
				mpfr_sub(y, y, MPFR_PI, MPFR_RNDN);
				mpfr_taylor_cos_bounded(R, y, N);
				mpfr_neg(R, R, MPFR_RNDN);
			}
		}
		else
		{
			if(mpfr_cmp(y, MPFR_HALF_PI) >= 0)
			{
				mpfr_sub(y, MPFR_PI, y, MPFR_RNDN);
				mpfr_taylor_cos_bounded(R, y, N);
				mpfr_neg(R, R, MPFR_RNDN);
			}
			else
			{
				mpfr_taylor_cos_bounded(R, y, N);
			}
		}
	}
	else
	{
		mpfr_neg(y, y, MPFR_RNDN);
		mpfr_taylor_cos_bounded(R, y, N);
	}
}
			
double taylor_sin(double x, unsigned int N)
{
	if(x >= 0)
	{
		while(x >= TWO_PI)
			x -= TWO_PI;

		if(x >= PI)
			if(x - PI >= HALF_PI)
				return -1 * taylor_sin_bounded(TWO_PI - x, N);
			else
				return -1 * taylor_sin_bounded(x - PI, N);
		else
			if(x >= HALF_PI)
				return taylor_sin_bounded(PI - x, N);
			else
				return taylor_sin_bounded(x, N);
	}

	return -1 * taylor_sin(-x, N);
}

void mpfr_taylor_sin(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_t y, t;
	mpfr_init_set(y, x, MPFR_RNDN);
	mpfr_init(t);

	if(mpfr_cmp_ui(y, 0) >= 0)
	{
		while(mpfr_cmp(y, MPFR_TWO_PI) >= 0)
			mpfr_sub(y, y, MPFR_TWO_PI, MPFR_RNDN);

		if(mpfr_cmp(y, MPFR_PI) >= 0)
		{
			mpfr_sub(t, y, MPFR_PI, MPFR_RNDN);
			if(mpfr_cmp(t, MPFR_PI) >= 0)
			{
				mpfr_sub(y, MPFR_TWO_PI, y, MPFR_RNDN);
				mpfr_taylor_sin_bounded(R, y, N);
				mpfr_neg(R, R, MPFR_RNDN);
			}
			else
			{
				mpfr_sub(y, y, MPFR_PI, MPFR_RNDN);
				mpfr_taylor_sin_bounded(R, y, N);
				mpfr_neg(R, R, MPFR_RNDN);
			}
		}
		else
		{
			if(mpfr_cmp(y, MPFR_HALF_PI) >= 0)
			{
				mpfr_sub(y, MPFR_PI, y, MPFR_RNDN);
				mpfr_taylor_sin_bounded(R, y, N);
			}
			else
			{
				mpfr_taylor_sin_bounded(R, y, N);
			}
		}
	}
	else
	{
		mpfr_neg(y, y, MPFR_RNDN);
		mpfr_taylor_sin(R, y, N);
		mpfr_neg(R, R, MPFR_RNDN);
	}
}
	

double taylor_tan(double x, unsigned int N)
{
	return taylor_sin(x,N)/taylor_cos(x,N);
}

void mpfr_taylor_tan(mpfr_t R, mpfr_t x, unsigned int N)
{
	mpfr_t S, C;
	mpfr_init(S);
	mpfr_init(C);
	mpfr_taylor_sin(S, x, N);
	mpfr_taylor_cos(C, x, N);
	mpfr_div(R, S, C, N);
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
					printf("Cos(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_cos(x, n));
				else
					printf("Usage: %s a <x=value for Cos(x)> <n> "
						   "<D=Number of digits to display>\n",

							argv[0]);
				break;
			 
			case 'b':
				if(argc == 5 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Sin(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_sin(x, n));
				else
					printf("Usage: %s b <x=value for Sin(x)> <n> "
						   "<D=Number of digits to display>\n",

							argv[0]);
				break;

			case 'c':
				if(argc == 5 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Tan(%.*lf) = %.*lf\n",
							d(D), x, D, taylor_tan(x, n));
				else
					printf("Usage: %s a <x=value for Tan(x)> <n> "
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

						sprintf(sf, "cos(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_cos(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Cos(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Cos(x)> "
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

						sprintf(sf, "sin(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_sin(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Sin(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Sin(x)> "
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

						sprintf(sf, "tan(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);
						
						mpfr_taylor_tan(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Tan(x)> "
							   "<D=Number of digits to display> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Tan(x)> "
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
