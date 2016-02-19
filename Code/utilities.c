#include <gmp.h>
#include <mpfr.h>
#include "utilities.h"

const double ROOT_2     = 1.4142135623730950488016887242096980785696718753;
const double ROOT_2_INV = 0.7071067811865475244008443621048490392848359376;

inline unsigned int d(unsigned int D)
{
	return D > 10 ? 10 : D;
}

void inline mpfr_digits_to_tolerance(unsigned int D, mpfr_t T)
{
	mpfr_init_set_ui(T, 10, MPFR_RNDN);
	mpfr_pow_ui(T, T, D, MPFR_RNDN);
	mpfr_ui_div(T, 1, T, MPFR_RNDN);
}
