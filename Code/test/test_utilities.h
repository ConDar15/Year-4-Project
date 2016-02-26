#ifndef TEST_UTILITIES_HEADER
#define TEST_UTILITIES_HEADER
	#include <time.h>
	
	clock_t time_average(clock_t*, int);
	clock_t time_total(clock_t*, int);
	clock_t time_min(clock_t*, int);
	clock_t time_max(clock_t*, int);
#endif
