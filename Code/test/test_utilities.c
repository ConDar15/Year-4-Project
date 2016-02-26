#include <time.h>

#include "test_utilities.h"

clock_t time_average(clock_t *times, int N)
{
	clock_t t = time_total(times, N);
	return t/N;
}

clock_t time_total(clock_t *times, int N)
{
	clock_t t = 0;
	for(int i = 0; i < N; i++)
		t += times[i];
	return t;
}

clock_t time_min(clock_t *times, int N)
{
	clock_t t = 0x0FFFFFFF;
	for(int i = 0; i < N; i++)
		if(times[i] < t)
			t = times[i];
	return t;
}

clock_t time_max(clock_t *times, int N)
{
	clock_t t = 0;
	for(int i = 0; i < N; i++)
		if(times[i] > t)
			t = times[i];
	return t;
}
