#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "geometric_trig.h"
#include "taylor_trig.h"
#include "cordic_trig.h"
#include "test_utilities.h"

int main(int argc, char **argv)
{
	FILE *in;
	int first = 1, c = 0;
	uintmax_t msec;
	unsigned int K, I;
	size_t r;
	char line[1024];
	clock_t *time[4], t;
	char *name[4] = {"geometric_cos", "taylor_cos", 
					 "cordic_cos", "builtin_cos"};

	double *N_d;
	//Assumes that a single argument is given and that it is a file of
	//	test values.
	sscanf(argv[2], "%u", &I);
	in = fopen(argv[1], "r");
	if (in == NULL)
	{
		printf("File %s not found for reading\n", argv[1]);
		exit(1);
	}
		
	while(fgets(line, sizeof(line), in))
	{
		if(first)
		{
			if(sscanf(line, "%u", &K) != 1)
			{
				printf("File %s incorrectly formatted\n", argv[1]);
				exit(2);
			}
			N_d = malloc(K * sizeof(*N_d));
			time[0] = malloc(K * sizeof(**time));
			time[1] = malloc(K * sizeof(**time));
			time[2] = malloc(K * sizeof(**time));
			time[3] = malloc(K * sizeof(**time));
			first = 0;
		}
		else
		{
			if(sscanf(line, "%lf", &N_d[c]) != 1)
				printf("Line %d is incorrectly formatted\n", c+2);
			else
				++c;
		}
		K = c;	
	} 
	fclose(in);
	int j;
	for(int i = 0; i < K; i++)
	{
		t = clock();
		for(j=0; j < I; j++)
			geometric_cos_bounded(N_d[i], 16);
		time[0][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			taylor_cos_bounded(N_d[i], 8);
		time[1][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			free(cordic_trig(N_d[i], 34));
		time[2][i] = clock() - t;

		t = clock();
		for(j=0; j < I; j++)
			cos(N_d[i]);
		time[3][i] = clock() - t;
	}
	printf("K = %u\n", K);
	for(int i = 0; i < 4; i++)
	{
		msec = (uintmax_t)time_average(time[i], K) * 1000 / CLOCKS_PER_SEC;
		printf("Average time for %s: %lu.%03lus\n",
		 	   name[i], msec/1000, msec % 1000);
		msec = (uintmax_t)time_total(time[i], K) * 1000 / CLOCKS_PER_SEC;
		printf("Total time for %s: %lu.%03lus\n",
			   name[i], msec/1000, msec % 1000);
		msec = (uintmax_t)time_min(time[i], K) * 1000 / CLOCKS_PER_SEC;
		printf("Minimum time for %s: %lu.%03lus\n",
			   name[i], msec/1000, msec % 1000);
		msec = (uintmax_t)time_max(time[i], K) * 1000 / CLOCKS_PER_SEC;
		printf("Maximum time for %s: %lu.%03lus\n",
			   name[i], msec/1000, msec % 1000);
		printf("\n");
	}

}
