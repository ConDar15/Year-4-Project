#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <assert.h>
#include <math.h>

#include "utilities.h"
#include "newton_inv_sqrt.h"

#define INIT_CONSTANTS mpfr_init_set_d(MPFR_THREE_HALF, 1.5, MPFR_RNDN); \
					   in = fopen(ROOT_2_INFILE, "r"); \
					   mpfr_init(MPFR_ROOT_2); \
					   mpfr_inp_str(MPFR_ROOT_2, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   in = fopen(ROOT_2_INV_INFILE, "r"); \
					   mpfr_init(MPFR_ROOT_2_INV); \
					   mpfr_inp_str(MPFR_ROOT_2_INV, in, 10, MPFR_RNDN); \
					   fclose(in);

mpfr_t MPFR_ROOT_2, MPFR_ROOT_2_INV, MPFR_THREE_HALF;

double newton_inv_sqrt(double N, double T)
{
	assert(N >= 0);
	assert(T >= 0);

	int e;
	double x, px, d, NN, N_2;

	NN  = N;
	N   = frexp(N, &e);
	N_2 = 0.5*N;

	x = 1;
	d = 1000000;

	while(d > T)
	{
		px = x;
		x = x * (1.5 - N_2*x*x);
		d = fabs(x - px);
	}

	if(e%2)
		x *= e > 0 ? ROOT_2_INV : ROOT_2;
	x *= NN;
	return ldexp(x, -e / 2);
}

void mpfr_newton_inv_sqrt(mpfr_t R, mpfr_t N, mpfr_t T)
{
	mpfr_t x, px, d, t, n, n_2;
	mpfr_exp_t e;

	mpfr_init(n);
	mpfr_frexp(&e, n, N, MPFR_RNDN);
	mpfr_init_set(n_2, n, MPFR_RNDN);
	mpfr_div_ui(n_2, n, 2, MPFR_RNDN);

	mpfr_init_set_ui(x, 1, MPFR_RNDN);
	mpfr_init(px);
	mpfr_init_set_ui(d, 1000000, MPFR_RNDN);
	mpfr_init(t);

	while(mpfr_cmp(d, T) > 0)
	{
		mpfr_set(px, x, MPFR_RNDN);
		mpfr_mul(t, x, x, MPFR_RNDN);
		mpfr_mul(t, t, n_2, MPFR_RNDN);
		mpfr_sub(t, MPFR_THREE_HALF, t, MPFR_RNDN);
		mpfr_mul(x, x, t, MPFR_RNDN);
		mpfr_sub(d, x, px, MPFR_RNDN);
		mpfr_abs(d, d, MPFR_RNDN);
	}

	if(e%2)
		if(e > 0)
			mpfr_mul(x, MPFR_ROOT_2_INV, x, MPFR_RNDN);
		else
			mpfr_mul(x, MPFR_ROOT_2, x, MPFR_RNDN);
	mpfr_mul(x, x, N, MPFR_RNDN);
	mpfr_mul_2si(R, x, -e/2, MPFR_RNDN);
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

	switch(argv[1][0])
	{
		case 'a':
			if (argc == 5 &&
					sscanf(argv[2], "%lf", &N) == 1 &&
					sscanf(argv[3], "%lf", &T) == 1 &&
					sscanf(argv[4], "%u" , &D) == 1)
				printf("sqrt(%.*lf) =~ %.*lf\n", d(D), N, D, newton_inv_sqrt(N, T));
			else
				printf("Usage: %s a <N=Value to sqrt> "
					   "<T=Tolerance> <D=Number of digits to display>\n",
					   argv[0]);
			break;

		case 'b':
			if (argc == 5 &&
					sscanf(argv[3], "%u", &D) == 1 &&
					sscanf(argv[4], "%u", &p) == 1)
			{
				mpfr_set_default_prec(p);
				INIT_CONSTANTS

				if(mpfr_init_set_str(Nr, argv[2], 10, MPFR_RNDN) == 0)
				{
					mpfr_init(R);

					mpfr_digits_to_tolerance(D, Tr);

					sprintf(sf, "sqrt(%%.%uRNf) =~\n\t%%.%uRNf\n", d(D), D);

					mpfr_newton_inv_sqrt(R, Nr, Tr);
					mpfr_printf(sf, Nr, R);
				}
				else
					printf("Usage: %s b <N=Value to sqrt> "
						   "<D=Number of digits to calculate to> "
						   "<p=bits of precision>\n", argv[0]);
			}
			else
				printf("Usage: %s b <N=Value to sqrt> "
					   "<D=Number of digits to calculate to> "
					   "<p=bits of precision>\n", argv[0]);
			break;

		default:
			printf("Usage: %s [a/b] <Arguments>\n", argv[0]);
	}
}
#endif
