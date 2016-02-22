#include <assert.h>
#include "trig_utilities.h"

TRIG_FIXED_TYPE double_to_fixed(double d)
{
	assert(d < 4 && d > -4);
	return (TRIG_FIXED_TYPE) (d * (TRIG_FIXED_TYPE)CONVERSION_VALUE);
}

double fixed_to_double(TRIG_FIXED_TYPE t)
{
	return (double)t / (TRIG_FIXED_TYPE)CONVERSION_VALUE;
}
