#ifndef HYPERBOLIC_LOG_HEADER
#define HYPERBOLIC_LOG_HEADER

	double hyperbolic_nat_log(double, unsigned int);
	double hyperbolic_log(double, double, unsigned int);
	void mpfr_hyperbolic_nat_log(mpfr_t, mpfr_t, mpz_t);
	void mpfr_hyperbolic_log(mpfr_t, mpfr_t, mpfr_t, mpz_t);

#endif
