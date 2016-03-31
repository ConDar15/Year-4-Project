#ifndef CONT_FRAC_EXP_HEADER
#define CONT_FRAC_EXP_HEADER

	double cont_frac_exp_v1(double, unsigned int);
	double cont_frac_exp_v2(double, unsigned int);
	double cont_frac_exp_v3(double, unsigned int);
	void mpfr_cont_frac_v3(mpfr_t, mpfr_t, mpz_t);

	double improved_pow(double, double, unsigned int);
	void mpfr_impfroved_pow(mpfr_t, mpfr_t, mpfr_t, mpz_t);

#endif
