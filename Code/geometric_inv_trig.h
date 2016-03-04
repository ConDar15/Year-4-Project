#ifndef GEOMETRIC_INV_TRIG_HEADER
#define GEOMETRIC_INV_TRIG_HEADER

	double geometric_acos_bounded(double, unsigned int);
	void mpfr_geometric_acos_bounded(mpfr_t, mpfr_t, unsigned int);
	
	double geometric_acos(double, unsigned int);
	double geometric_asin(double, unsigned int);
	double geometric_atan(double, unsigned int);

	void mpfr_geometric_acos(mpfr_t, mpfr_t, unsigned int);
	void mpfr_geometric_asin(mpfr_t, mpfr_t, unsigned int);
	void mpfr_geometric_atan(mpfr_t, mpfr_t, unsigned int);
	
#endif
