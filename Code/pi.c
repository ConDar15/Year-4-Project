#include <gmp.h>
#include <mpfr.h>
#include <stdio.h>

void mpfr_pi_Chudnovsky(mpfr_t pi, unsigned int N)
{
	mpz_t a, b, d, e, f;
	mpfr_t c, t;

	mpz_inits(a, b, d, e, f, NULL);
	mpfr_init(t);
	mpfr_init(c);
	
	mpfr_set_ui(pi, 0, MPFR_RNDN);
	for(unsigned int k = 0; k < N; k++)
	{
		mpz_fac_ui(a, 3*k);
		
		mpz_fac_ui(b, k);
		mpz_pow_ui(b, b, 3);
		
		mpfr_set_ui(c, 3*k, MPFR_RNDN);
		mpfr_add_d(c, c, 1.5, MPFR_RNDN);
		mpfr_ui_pow(c, 640320, c, MPFR_RNDN);
		
		mpz_fac_ui(d, 6*k);
		
		mpz_set_ui(e, 545140134);
		mpz_mul_ui(e, e, k);
		mpz_add_ui(e, e, 13591409);
		
		mpz_set_si(f, -1);
		mpz_pow_ui(f, f, k);

		mpz_mul(a, a, b);
		mpz_mul(d, d, e);
		mpz_mul(d, d, f);

		mpfr_set_z(t, d, MPFR_RNDN);
		mpfr_div_z(t, t, a, MPFR_RNDN);
		mpfr_div(t, t, c, MPFR_RNDN);
		mpfr_add(pi, pi, t, MPFR_RNDN);
	}
	mpfr_mul_ui(pi, pi, 12, MPFR_RNDN);
	mpfr_ui_div(pi, 1, pi, MPFR_RNDN);
}

int main(int argc, char **argv)
{
	mpfr_t pi;
	unsigned int N, P, D;
	char sf[50];

	sscanf(argv[1], "%u", &N);
	sscanf(argv[2], "%u", &D);
	sscanf(argv[3], "%u", &P);
	mpfr_set_default_prec(P);
	mpfr_init(pi);
	mpfr_pi_Chudnovsky(pi, N);
	sprintf(sf, "%%.%uRZf", D);
	mpfr_printf(sf, pi);
}
