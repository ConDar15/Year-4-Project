#ifndef TAYLOR_TRIG_HEADER
#define TAYLOR_TRIG_HEADER

	double taylor_cos_bounded(double, unsigned int);
	double taylor_sin_bounded(double, unsigned int);
	double taylor_sin(double, unsigned int);
	double taylor_cos(double, unsigned int);
	double taulor_tan(double, unsigned int);

	void mpfr_taylor_cos_bounded(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_sin_bounded(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_sin(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_cos(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_cos(mpfr_t, mpfr_t, unsigned int);

#endif
