#ifndef TAYLOR_INV_TRIG_HEADER
#define TAYLOR_INV_TRIG_HEADER
	
	double taylor_asin(double, unsigned int);
	double taylor_acos(double, unsigned int);
	double taylor_atan_bounded(double, unsigned int);
	double taylor_atan(double, unsigned int);

	void mpfr_taylor_asin(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_acos(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_atan_bounded(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_atan(mpfr_t, mpfr_t, unsigned int);
#endif
