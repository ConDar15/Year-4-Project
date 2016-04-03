#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

#include "log_exp_utilities.h"
#include "utilities.h"
#include "hyperbolic_log.h"
#include "int_exp.h"
#include "cont_frac_exp.h"

#define INIT_CONSTANTS in = fopen(NAT_LOG_2_INFILE, "r"); \
					   mpfr_init(MPFR_NAT_LOG_2); \
					   mpfr_inp_str(MPFR_NAT_LOG_2, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   in = fopen(E_CONST_INFILE, "r"); \
					   mpfr_init(MPFR_E_CONST); \
					   mpfr_inp_str(MPFR_E_CONST, in, 10, MPFR_RNDN); \
					   fclose(in); \
					   mpfr_init_set_ui(MPFR_TWO, 2, MPFR_RNDN);

mpfr_t MPFR_NAT_LOG_2, MPFR_TWO, MPFR_E_CONST;

double cont_frac_exp_v1(double x, unsigned int n)
{
	double a = x, pA, A, nA, ca, pB, B, nB, cb;
	unsigned int b = 0;

	if (x < 0)
		//Calculates the reciprocal if x < 0
		return 1/cont_frac_exp_v1(-x, n);
	else if (x == 0)
		//Basic value check
		return 1;
	else if (x > 1)
		//1/2 <= a < 1 and x == a*2^b
		a = frexp(x, &b);
	
	//Initialises the previous values
	//	pA = A_1, pB = B_1
	pA = a + 1;
	pB = 1;
	//Initialises the current values
	//	A = A_2, B = B_2
	A  = a*a + 2*a + 2;
	B  = 2;
	
	//Initialises the co-efficients
	ca = -a;
	cb = 2 + a;
	
	for(unsigned int k = 2; k <= n; ++k)
	{
		ca -= a;
		++cb;
		
		//Calculates the updates to A and B
		nA = cb*A + ca*pA;
		nB = cb*B + ca*pB;
		
		//Ensures the variables hold the correct values
		pA = A;
		pB = B;
		A  = nA;
		B  = nB;
	}
	
	//If b b==0 then return A/B, otherwise return (A/B)^(2^b)
	return b ? squaring_int_exp(A/B, (unsigned int)squaring_int_exp(2, b))
			 : A/B;
}

double cont_frac_exp_v2(double x, unsigned int n)
{
	double a = x, pA, A, nA, ca, pB, B, nB;
	unsigned int b = 0, cb;
	
	if (x < 0)
		//Calculates the reciprocal if x < 0
		return 1/cont_frac_exp_v1(-x, n);
	else if (x == 0)
		//Basic value check
		return 1;
	else if (x > 1)
		//1/2 <= a < 1 and x == a*2^b
		a = frexp(x, &b);

	//Sets the initial values
	//	pA = A_1, pB = B_1
	//	A  = A_2,  B = B_2	
	pA = 1;
	pB = 1;
	A  = 1;
	B  = 1 - a;

	//Initialises the co-efficients
	ca = 0;
	cb = 1;

	//Main Loop
	for(unsigned int k = 3; k <= n; ++k)
	{
		++cb;
	
		if (k % 2)
		{
			ca += a;
			nA = cb*A + ca*pA;
			nB = cb*B + ca*pB;
		}
		else
		{
			nA = cb*A - ca*pA;
			nB = cb*B - ca*pB;
		}
			
		pA = A;
		pB = B;
		A  = nA;
		B  = nB;
	}

	//If b b==0 then return A/B, otherwise return (A/B)^(2^b)
	return b ? squaring_int_exp(A/B, (unsigned int)squaring_int_exp(2, b))
			 : A/B;
}

double cont_frac_exp_v3(double x, unsigned int n)
{
	double a = x, pA, A, nA, ca, pB, B, nB;
	unsigned int b = 0, cb;
	
	if (x < 0)
		//Calculates the reciprocal if x < 0
		return 1/cont_frac_exp_v1(-x, n);
	else if (x == 0)
		//Basic value check
		return 1;
	else if (x > 1)
		//1/2 <= a < 1 and x == a*2^b
		a = frexp(x, &b);
		
	//Sets the initial values
	//	pA = A_0, pB = B_0
	//	A  = A_1,  B = B_1	
	pA = 1;
	pB = 1;
	A  = 2 + a;
	B  = 2 - a;

	//Initialises the co-efficients
	ca = a*a;
	cb = 2;

	//Main loop
	for(unsigned int k = 2; k <= n; ++k)
	{
		cb += 4;
		
		nA = cb*A + ca*pA;
		nB = cb*B + ca*pB;
		
		pA = A;
		pB = B;
		A  = nA;
		B  = nB;
	}

	//If b b==0 then return A/B, otherwise return (A/B)^(2^b)
	return b ? squaring_int_exp(A/B, (unsigned int)squaring_int_exp(2, b))
			 : A/B;
}

void mpfr_cont_frac_exp_v3(mpfr_t R, mpfr_t x, mpz_t n)
{
	assert(mpz_cmp_ui(n, 0) >= 0);

	mpfr_t a, pA, A, nA, ca, pB, B, nB, C, tmp1, tmp2;
	mpfr_exp_t b = 0;
	mpz_t N, k, cb;

	if (mpfr_cmp_ui(x, 0) < 0)
	{
		mpfr_neg(x, x, MPFR_RNDN);
		mpfr_cont_frac_exp_v3(R, x, n);
		mpfr_ui_div(R, 1, R, MPFR_RNDN);
	}
	else if (mpfr_cmp_ui(x, 0) == 0)
	{
		mpfr_set_ui(R, 1, MPFR_RNDN);
	}
	else
	{
		mpfr_init(a);
		if (mpfr_cmp_ui(x, 1) > 0)
		{
			mpfr_frexp(&b, a, x, MPFR_RNDN);
			mpz_init(N);
			mpfr_init(C);
			mpz_ui_pow_ui(N, 2, b);
		}
		else
		{
			mpfr_init_set_ui(C, 0, MPFR_RNDN);
			mpfr_set(a, x, MPFR_RNDN);
		}
		
		mpfr_init_set_ui(pA, 1, MPFR_RNDN);
		mpfr_init_set_ui(pB, 1, MPFR_RNDN);
		mpfr_init_set_ui(A,  2, MPFR_RNDN);
		mpfr_init_set_ui(B,  2, MPFR_RNDN);
		mpfr_add(A, A, a, MPFR_RNDN);
		mpfr_sub(B, B, a, MPFR_RNDN);

		mpfr_init_set(ca, a, MPFR_RNDN);
		mpfr_mul(ca, ca, a, MPFR_RNDN);
		mpz_init_set_ui(cb, 2);

		mpfr_init(nA);
		mpfr_init(nB);
		mpfr_init(tmp1);
		mpfr_init(tmp2);

		for(mpz_init_set_ui(k, 2); mpz_cmp(k, n) <= 0; mpz_add_ui(k, k, 1))
		{
			mpz_add_ui(cb, cb, 4);

			mpfr_mul_z(tmp1, A, cb, MPFR_RNDN);
			mpfr_mul(tmp2, ca, pA, MPFR_RNDN);
			mpfr_add(nA, tmp1, tmp2, MPFR_RNDN);

			mpfr_mul_z(tmp1, B, cb, MPFR_RNDN);
			mpfr_mul(tmp2, ca, pB, MPFR_RNDN);
			mpfr_add(nB, tmp1, tmp2, MPFR_RNDN);

			mpfr_set(pA, A, MPFR_RNDN);
			mpfr_set(pB, B, MPFR_RNDN);
			mpfr_set(A, nA, MPFR_RNDN);
			mpfr_set(B, nB, MPFR_RNDN);
		}

		mpfr_div(R, A, B, MPFR_RNDN);
		
		if(b)
		{
			mpfr_set(A, R, MPFR_RNDN);
			mpfr_squaring_int_exp(R, A, N);
		}
	}
}

double improved_pow(double x, double y, unsigned int n)
{
	assert(x > 0);
	
	return cont_frac_exp_v3(y * hyperbolic_nat_log(x, n), n);
}

void mpfr_improved_pow(mpfr_t R, mpfr_t x, mpfr_t y, mpz_t n)
{
	assert(mpz_cmp_ui(n, 0) >= 0);
	assert(mpfr_cmp_ui(x, 0) > 0);

	mpfr_t A;
	mpfr_init(A);

	mpfr_hyperbolic_nat_log(A, x, n);
	mpfr_mul(A, y, A, MPFR_RNDN);
	mpfr_cont_frac_exp_v3(R, A, n);
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
					printf("exp(%.*lf) ~= %.*lf (v1)\n",
						d(D), x, D, cont_frac_exp_v1(x, n));
				else
					printf("Usage: %s a <x=Value for exp(x)> <n> "
						   "<D=Digits to display\n", argv[0]);
				break;
			
			case 'b':
				if(argc == 5 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%u",  &n) == 1 &&
					sscanf(argv[4], "%u",  &D) == 1)
					printf("exp(%.*lf) ~= %.*lf (v2)\n",
						d(D), x, D, cont_frac_exp_v2(x, n));
				else
					printf("Usage: %s b <x=Value for exp(x)> <n> "
						   "<D=Digits to display\n", argv[0]);
				break;

			case 'c':
				if(argc == 5 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%u",  &n) == 1 &&
					sscanf(argv[4], "%u",  &D) == 1)
					printf("exp(%.*lf) ~= %.*lf (v3)\n",
						d(D), x, D, cont_frac_exp_v3(x, n));
				else
					printf("Usage: %s c <x=Value for exp(x)> <n> "
						   "<D=Digits to display\n", argv[0]);
				break;

			case 'd':
				if(argc == 6 &&
					sscanf(argv[2], "%lf", &x) == 1 &&
					sscanf(argv[3], "%lf", &y) == 1 &&
					sscanf(argv[4], "%u",  &n) == 1 &&
					sscanf(argv[5], "%u",  &D) == 1)
					printf("pow(%.*lf, %.*lf) ~= %.*lf\n",
						d(D), x, d(D), y, D, improved_pow(x, y, n));
				else
					printf("Usage: %s d <x=Value for pow(x,y)> "
						   "<y=Value for pow(x,y)> <n> "
						   "<D=Digits to display>\n", argv[0]);
				break;

			case 'e':
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

						sprintf(sf, "exp(%%.%uRNf) =~\n\t%%.%uRNf\n",
							d(D), D);

						mpfr_cont_frac_exp_v3(R, X, N);
						mpfr_printf(sf, X, R);
					}
					else
						printf("Usage: %s e <X=Value for exp(X)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s e <X=Value for exp(X)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;

			case 'f':
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

						sprintf(sf, "pow(%%.%uRNf, %%.%uRNf) =~"
								    "\n\t%%.%uRNf\n",
								d(D), d(D), D);

						mpfr_improved_pow(R, X, Y, N);
						mpfr_printf(sf, X, Y, R);
					}
					else
						printf("Usage: %s f <X=Value for pow(X, Y)> "
							   "<Y=Value for pow(X,Y)> <N> "
							   "<D=Digits to display> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s f <X=Value for pow(X, Y)> "
						   "<Y=Value for pow(X,Y)> <N> "
						   "<D=Digits to display> "
						   "<p=bits of precision>\n", argv[0]);
				break;
			
			default:
				printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
}
#endif
