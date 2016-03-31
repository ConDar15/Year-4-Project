#ifndef UTILITIES_HEADER
	#define UTILITIES_HEADER

	#include <mpfr.h>

	#define ROOT_2_INFILE "root_2_digits.txt"
	#define ROOT_2_INV_INFILE "root_2_inv_digits.txt"

	extern const double ROOT_2, ROOT_2_INV;
	extern mpfr_t MPFR_ROOT_2, MPFR_ROOT_2, 
				  MPFR_ONE, MPFR_HALF, MPFR_THREE_HALF, MPFR_TWO; 

	unsigned int d(unsigned int);
	void mpfr_digits_to_tolerance(unsigned int, mpfr_t);
#endif
