#ifndef BISECT_ROOT_HEADER
	#define BISECT_ROOT_HEADER

	double bisect_sqrt(double, double);
	double bisect_sqrt_it(double, unsigned int);
	double ipow(double, unsigned int);
	double bisect_nRoot(double, double, unsigned int);
	void mpfr_bisect_sqrt(mpfr_t, mpfr_t, mpfr_t);
	void mpfr_bisect_nRoot(mpfr_t, mpfr_t, mpfr_t, unsigned int);
#endif
