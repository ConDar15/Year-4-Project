#ifndef INT_EXP_HEADER
#define INT_EXP_HEADER

	double naive_int_exp(const double, const int);
	double squaring_int_exp(const double, const int);
	void mpfr_naive_int_exp(mpfr_t, mpfr_t, mpz_t);
	void mpfr_squaring_int_exp(mpfr_t, mpfr_t, mpz_t);

#endif
