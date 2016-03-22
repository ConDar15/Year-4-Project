#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#include "int_exp.h"
#include "utilities.h"

double naive_int_exp(const double x, const int a)
{
	if(a < 0)
		return 1/naive_int_exp(x, -a);
	double z = 1;
	int n = a;
	while(n--)
		z *= x;
	return z;
}

double squaring_int_exp(const double x, const int a)
{
	if(a < 0)
		return 1/squaring_int_exp(x, -a);
	double y = x, z = 1;
	int n = a;
	while(n)
	{
		if(n%2)
		{
			z *= y;
			--n;
		}
		y *= y;
		n >>= 1;
	}
	return z;
}

void mpfr_naive_int_exp(mpfr_t R, mpfr_t x, mpz_t a)
{
	if(mpz_cmp_ui(a, 0) < 0)
	{
		mpz_t b;
		mpz_init_set(b, a);
		mpz_neg(b, b);
		mpfr_naive_int_exp(R, x, b);
		mpfr_ui_div(R, 1, R, MPFR_RNDN);
	}
	else
	{
		mpz_t n;
		mpz_init_set(n, a);
		mpfr_set_ui(R, 1, MPFR_RNDN);
		while(mpz_cmp_ui(n, 0) > 0)
		{
			mpfr_mul(R, R, x, MPFR_RNDN);
			mpz_sub_ui(n, n, 1);
		}
	}
}

void mpfr_squaring_int_exp(mpfr_t R, mpfr_t x, mpz_t a)
{
	if(mpz_cmp_ui(a, 0) < 0)
	{
		mpz_t b;
		mpz_init_set(b, a);
		mpz_neg(b, b);
		mpfr_squaring_int_exp(R, x, b);
		mpfr_ui_div(R, 1, R, MPFR_RNDN);
	}
	else
	{
		mpfr_t y;
		mpz_t n;
		mpfr_init_set(y, x, MPFR_RNDN);
		mpfr_set_ui(R, 1, MPFR_RNDN);
		mpz_init_set(n, a);
		while(mpz_cmp_ui(n, 0) > 0)
		{
			if(mpz_odd_p(n))
			{
				mpfr_mul(R, R, y, MPFR_RNDN);
				mpz_sub_ui(n, n, 1);
			}
			mpfr_mul(y, y, y, MPFR_RNDN);
			mpz_div_ui(n, n, 2);
		}
	}
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double x;
	unsigned int n, D, p;
	mpfr_t X, R;
	mpz_t N;
	char sf[50];

	if(argc > 1)
	{
		switch(argv[1][0])
		{
			case 'a':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("(%.*lf)^(%d) = %.*lf (Naive)\n",
							d(D), x, n, D, naive_int_exp(x, n));
				else
					printf("Uasge: %s a <x=Base for exp> "
					       "<n=Exponent for exp> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;

			case 'b':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("(%.*lf)^(%d) = %.*lf (Squaring)\n",
							d(D), x, n, D, squaring_int_exp(x, n));
				else
					printf("Uasge: %s b <x=Base for exp> "
					       "<n=Exponent for exp> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;

			case 'c':
				if (argc == 6 &&
						sscanf(argv[4], "%u", &D) == 1 &&
						sscanf(argv[5], "%u" , &p) == 1)
				{
					mpfr_set_default_prec(p);
								
					if (mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN)==0 &&
						mpz_init_set_str(N, argv[3], 10) == 0)
					{
						mpfr_init(R);
					
						sprintf(sf, "(%%.%uRNf)^%%Zd =~\t(Naive)"
									"\n\t%%.%uRNf\n", d(D), D);
					
						mpfr_naive_int_exp(R, X, N);
						mpfr_printf(sf, X, N, R);
					}
					else
						printf("Usage: %s c <X=Base for exp> "
						       "<N=Exponent for exp> "
							   "<D=Number of digitsto calculate to> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s c <X=Base for exp> "
					       "<N=Exponent for exp> "
						   "<D=Number of digitsto calculate to> "
						   "<p=bits of precision>\n", argv[0]);
				break;

			case 'd':
				if (argc == 6 &&
						sscanf(argv[4], "%u", &D) == 1 &&
						sscanf(argv[5], "%u" , &p) == 1)
				{
					mpfr_set_default_prec(p);
								
					if (mpfr_init_set_str(X, argv[2], 10, MPFR_RNDN)==0 &&
						mpz_init_set_str(N, argv[3], 10) == 0)
					{
						mpfr_init(R);
					
						sprintf(sf, "(%%.%uRNf)^%%Zd =~\t(Squaring)"
									"\n\t%%.%uRNf\n", d(D), D);
						
						mpfr_squaring_int_exp(R, X, N);
						mpfr_printf(sf, X, N, R);
					}
					else
						printf("Usage: %s d <X=Base for exp> "
						       "<N=Exponent for exp> "
							   "<D=Number of digitsto calculate to> "
							   "<p=bits of precision>\n", argv[0]);
				}
				else
					printf("Usage: %s d <X=Base for exp> "
					       "<N=Exponent for exp> "
						   "<D=Number of digitsto calculate to> "
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
