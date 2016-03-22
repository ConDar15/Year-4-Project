#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

#include "exact_root.h"
#include "bisect_root.h"
#include "newton_sqrt.h"
#include "newton_inv_sqrt.h"
#include "test_utilities.h"

int main(int argc, char **argv)
{
	FILE *in;
	int first = 1, c = 0, NUM = 5;
	uintmax_t msec;
	unsigned int K, I;
	size_t r;
	char line[1024];
	clock_t *time[5], t;
	char *name[] = {"root_digits_precise", "bisect_sqrt", 
					 "newton_sqrt", "newton_inv_sqrt", "builtin_sqrt"};

	double *N_d;
	char **L;
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
			L   = malloc(K * sizeof(*L));
			for(int i = 0; i < NUM; ++i)
				time[i] = malloc(K * sizeof(**time));
			for(int i = 0; i < K; ++i)
				L[i] = malloc(32 * sizeof(**L));
			first = 0;
		}
		else
		{
			if(sscanf(line, "%lf", &N_d[c]) != 1)
				printf("Line %d is incorrectly formatted\n", c+2);
			else
				strcpy(L[c++], line);
		}
		K = c;	
	} 
	fclose(in);
	int j;
	for(int i = 0; i < K; i++)
	{
		t = clock();
		for(j=0; j < I; j++)
			free(root_digits_precise(L[i], 11));
		time[0][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			bisect_sqrt_it(N_d[i], 33);
		time[1][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			newton_sqrt_v3_it(N_d[i], 5);
		time[2][i] = clock() - t;

		t = clock();
		for(j=0; j < I; j++)
			newton_inv_sqrt_it(N_d[i], 5);
		time[3][i] = clock() - t;

		t = clock();
		for(j=0; j < I; j++)
			sqrt(N_d[i]);
		time[4][i] = clock() - t;
	}
	for(int i = 0; i < NUM; i++)
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
