#ifndef GEOMETRIC_TRIG_HEADER
#define GEOMETRIC_TRIG_HEADER

	double geometric_cos_bounded(double, unsigned int);
	void mpfr_geometric_cos_bounded(mpfr_t, mpfr_t, unsigned int);
	
	double geometric_cos(double, unsigned int);
	double geometric_sin(double, unsigned int);
	double geometric_tan(double, unsigned int);

	void mpfr_geometric_cos(mpfr_t, mpfr_t, unsigned int);
	void mpfr_geometric_sin(mpfr_t, mpfr_t, unsigned int);
	void mpfr_geometric_tan(mpfr_t, mpfr_t, unsigned int);

#endif
