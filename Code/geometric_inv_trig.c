#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#include "geometric_inv_trig.h"
#include "trig_utilities.h"
#include "utilities.h"

#define INIT_CONSTANTS in = fopen(PI_INFILE, "r"); \
					   mpfr_init(MPFR_PI); \
					   mpfr_inp_str(MPFR_PI, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   mpfr_init(MPFR_HALF_PI); \
					   mpfr_div_ui(MPFR_HALF_PI, MPFR_PI, 2, MPFR_RNDN);

mpfr_t MPFR_PI, MPFR_HALF_PI;

double geometric_acos_bounded(double x, unsigned int n)
{
	//Ensures the given value is a valid cosine value
	assert(x >= 0 && x <= 1);
	
	//Reversing the last line of geometric_cos_bounded
	double h = 2-2*x;
	
	//Reverses the iterative process oof geometric_cos_bounded
	for(int i = 0; i < n; i++)
		h = 2 - sqrt(4 - h);

	//Reverses the initialisation proceduce in geometric_cos_bounded
	h *= pow(4, n);
	return sqrt(h);
}

void mpfr_geometric_acos_bounded(mpfr_t R, mpfr_t x, unsigned int n)
{
	assert(mpfr_cmp_ui(x, 0) >= 0 && mpfr_cmp_ui(x, 1) <= 0);

	mpfr_t h;

	mpfr_init(h);
	mpfr_ui_sub(h, 1, x, MPFR_RNDN);
	mpfr_mul_ui(h, h, 2, MPFR_RNDN);
	
	for(int i = 0; i < n; i++)
	{
		mpfr_ui_sub(h, 4, h, MPFR_RNDN);
		mpfr_sqrt(h, h, MPFR_RNDN);
		mpfr_ui_sub(h, 2, h, MPFR_RNDN);
	}
	
	mpfr_ui_pow_ui(R, 4, n, MPFR_RNDN);
	mpfr_mul(R, h, R, MPFR_RNDN);
	mpfr_sqrt(R, R, MPFR_RNDN);
}

double geometric_acos(double x, unsigned int n)
{
	assert(x >= -1 && x <= 1);
	return x >= 0 ? geometric_acos_bounded(x,n)
				  : PI - geometric_acos_bounded(-x,n);
}

void mpfr_geometric_acos(mpfr_t R, mpfr_t x, unsigned int n)
{
	assert(mpfr_cmp_si(x, -1) >= 0 && mpfr_cmp_si(x, 1) <= 0);
	mpfr_t y;

	if(mpfr_cmp_ui(x, 0) < 0)
	{
		mpfr_init_set(y, x, MPFR_RNDN);
		mpfr_neg(y, x, MPFR_RNDN);
		mpfr_geometric_acos_bounded(R,y,n);
		mpfr_sub(R, MPFR_PI, R, MPFR_RNDN);
	}
	else
		mpfr_geometric_acos_bounded(R,x,n);
}

double geometric_asin(double x, unsigned int n)
{
	assert(x >= -1 && x <= 1);
	return HALF_PI - geometric_acos(x, n);
}

void mpfr_geometric_asin(mpfr_t R, mpfr_t x, unsigned int n)
{
	assert(mpfr_cmp_si(x, -1) >= 0 && mpfr_cmp_si(x, 1) <= 0);
	
	mpfr_geometric_acos(R, x, n);
	mpfr_sub(R, MPFR_HALF_PI, R, MPFR_RNDN);
}

double geometric_atan(double x, unsigned int n)
{
	return geometric_asin(x/sqrt(x*x + 1), n);
}

void mpfr_geometric_atan(mpfr_t R, mpfr_t x, unsigned int n)
{
	mpfr_t y;
	
	mpfr_init(y);
	mpfr_mul(y, x, x, MPFR_RNDN);
	mpfr_add_ui(y, y, 1, MPFR_RNDN);
	mpfr_sqrt(y, y, MPFR_RNDN);
	mpfr_div(y, x, y, MPFR_RNDN);
	
	mpfr_geometric_asin(R, y, n);
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
							d(D), x, D, geometric_acos(x, n));
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
							d(D), x, D, geometric_asin(x, n));
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
							d(D), x, D, geometric_atan(x, n));
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
						
						mpfr_geometric_acos(R, X, n);
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
						
						mpfr_geometric_asin(R, X, n);
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
						
						mpfr_geometric_atan(R, X, n);
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

