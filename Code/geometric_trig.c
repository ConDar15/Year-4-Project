#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#include "geometric_trig.h"
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

double geometric_cos_bounded(double x, unsigned int n)
{
	//Ensures that x is in the range [0, HALF_PI) and raises an error
	//	message if this is not the case.
	assert(x >= 0 && x < HALF_PI);
	
	//Sets the first chord length that will be the basis or our induction
	double h = (x*x)/pow(4, n);
	
	//Performs the induction steps
	for(int i = 0; i < n; i++)
		h = h*(4-h);
	//Returns the approximation of cos(x)
	return 1 - h/2;
}

void mpfr_geometric_cos_bounded(mpfr_t R, mpfr_t x, unsigned int n)
{
	assert(mpfr_cmp_ui(x, 0) >= 0 && mpfr_cmp(x, MPFR_HALF_PI) < 0);

	mpfr_t h, t; 
	mpz_t k;

	mpfr_init(t);

	mpfr_init(h);
	mpfr_mul(h, x, x, MPFR_RNDN);
	mpz_init(k);
	mpz_ui_pow_ui(k, 4, n);
	mpfr_div_z(h, h, k, MPFR_RNDN);


	for(int i = 0; i < n; i++)
	{
		mpfr_ui_sub(t, 4, h, MPFR_RNDN);
		mpfr_mul(h, h, t, MPFR_RNDN);
	}

	mpfr_div_ui(h, h, 2, MPFR_RNDN);
	mpfr_ui_sub(R, 1, h, MPFR_RNDN);
}


double geometric_cos(double x, unsigned int n)
{
	// We have two cases to consider, x >= 0 and x < 0
	if(x >= 0)
	{
		//Ensures x is in the range [0, TWO_PI) as ;
		//	cos(x + TWO_PI) == cos(x)
		while(x >= TWO_PI)
			x -= TWO_PI;

		//Calcualtes the correct modification of x to accurately
		//	calculate cos(x) when it is reduced to the range [0, HALF_PI)
		if(x >= PI)
			if(x - PI >= HALF_PI)
				return geometric_cos_bounded(TWO_PI - x, n);
			else
				return -1 * geometric_cos_bounded(x - PI, n);		
		else
			if(x >= HALF_PI)
				return -1 * geometric_cos_bounded(PI - x, n);
			else
				return geometric_cos_bounded(x, n);
	}
	
	//The second case increases x until it is non-negative and then
	//	recursively calls geometric_cos, to perform the code above
	while(x < 0)
		x += TWO_PI;
	return geometric_cos(x, n);
}

void mpfr_geometric_cos(mpfr_t R, mpfr_t x, unsigned int n)
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
				mpfr_geometric_cos_bounded(R, y ,n);
			}
			else
			{
				mpfr_sub(y, y, MPFR_PI, MPFR_RNDN);
				mpfr_geometric_cos_bounded(R, y, n);
				mpfr_mul_si(R, R, -1, MPFR_RNDN);
			}
		}
		else
		{	
			if(mpfr_cmp(y, MPFR_HALF_PI) >= 0)
			{
				mpfr_sub(y, MPFR_PI, y, MPFR_RNDN);
				mpfr_geometric_cos_bounded(R, y, n);
				mpfr_mul_si(R, R, -1, MPFR_RNDN);
			}
			else
			{
				mpfr_geometric_cos_bounded(R, y, n);
			}
		}
	}
	else
	{
		while(mpfr_cmp_ui(y, 0) < 0)
			mpfr_add(y, y, MPFR_TWO_PI, MPFR_RNDN);
		mpfr_geometric_cos(R, y, n);
	}
}

//sin(x) = cos(x - HALF_PI)
double geometric_sin(double x, unsigned int n)
{
	return geometric_cos(x - HALF_PI, n);
}

void mpfr_geometric_sin(mpfr_t R, mpfr_t x, unsigned int n)
{
	mpfr_t y;
	mpfr_init(y);
	mpfr_sub(y, x, MPFR_HALF_PI, MPFR_RNDN);
	mpfr_geometric_cos(R, y, n);
}

//tan(x) = sin(x)/cos(x)
double geometric_tan(double x, unsigned int n)
{
	return geometric_sin(x, n)/geometric_cos(x, n);
}

void mpfr_geometric_tan(mpfr_t R, mpfr_t x, unsigned int n)
{
	mpfr_t S, C;
	mpfr_init(S);
	mpfr_init(C);
	mpfr_geometric_sin(S, x, n);
	mpfr_geometric_cos(C, x, n);
	mpfr_div(R, S, C, MPFR_RNDN);
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
				if(argc == 4 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1)
					printf("Cos(%.*lf) = %.*lf\n",
							n+5, x, n+5, geometric_cos(x, n));
				else
					printf("Usage: %s a <x=value for cos(x)>, <n>\n",
							argv[0]);
				break;
			 
			case 'b':
				if(argc == 4 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1)
					printf("Sin(%.*lf) = %.*lf\n",
							n+5, x, n+5, geometric_sin(x, n));
				else
					printf("Usage: %s b <x=value for sin(x)>, <n>\n",
							argv[0]);
				break;

			case 'c':
				if(argc == 4 && 
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1)
					printf("Tan(%.*lf) = %.*lf\n",
							n+5, x, n+5, geometric_tan(x, n));
				else
					printf("Usage: %s a <x=value for tan(x)>, <n>\n",
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
						
						mpfr_geometric_cos(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Cos(x)> "
							   "<D=Number of digits to display<n> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Cos(x)> "
						   "<D=Number of digits to display<n> "
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
						
						mpfr_geometric_sin(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Sin(x)> "
							   "<D=Number of digits to display<n> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Sin(x)> "
						   "<D=Number of digits to display<n> "
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
						
						mpfr_geometric_tan(R, X, n);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s d <x=value for Tan(x)> "
							   "<D=Number of digits to display<n> "
							   "<n> <p=bits of precision to use>\n",
							   argv[0]);
				}
				else
					printf("Usage: %s d <x=value for Tan(x)> "
						   "<D=Number of digits to display<n> "
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
