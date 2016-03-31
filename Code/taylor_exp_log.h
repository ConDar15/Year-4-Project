#ifndef TAYLOR_EXP_LOG_HEADER
#define TAYLOR_EXP_LOG_HEARDR

	double naive_exp(double, unsigned int);
	void mpfr_naive_exp(mpfr_t, mpfr_t, mpz_t);
	
	double taylor_exp(double, unsigned int);
	double taylor_nat_log(double, unsigned int);

	double taylor_log(double, double, unsigned int);
	double taylor_pow(double, double, unsigned int);

	void mpfr_taylor_exp(mpfr_t, mpfr_t, mpz_t);
	void mpfr_taylor_nat_log(mpfr_t, mpfr_t, mpz_t);

	void mpfr_taylor_log(mpfr_t, mpfr_t, mpfr_t, mpz_t);
	void mpfr_taylor_pow(mpfr_t, mpfr_t, mpfr_t, mpz_t);

#endif
