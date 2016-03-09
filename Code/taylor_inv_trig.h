#ifndef TAYLOR_INV_TRIG_HEADER
#define TAYLOR_INV_TRIG_HEADER
	
	double taylor_aSin(double, unsigned int);
	double taylor_aCos(double, unsigned int);
	double taylor_aTan_bounded(double, unsigned int);
	double taylor_aTan(double, unsigned int);

	void mpfr_taylor_aSin(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_aCos(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_aTan_bounded(mpfr_t, mpfr_t, unsigned int);
	void mpfr_taylor_aTan(mpfr_t, mpfr_t, unsigned int);
#endif
