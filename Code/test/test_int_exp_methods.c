#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>

#include "int_exp.h"
#include "test_utilities.h"

int main(int argc, char **argv)
{
	FILE *in;
	int first = 1, c = 0, NUM = 2;
	uintmax_t msec;
	unsigned int K, I, *A_u;
	size_t r;
	char line[1024];
	clock_t *time[2], t;
	char *name[] = {"naive_int_exp", "squaring_int_exp"};

	double *X_d;
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
			X_d = malloc(K * sizeof(*X_d));
			A_u = malloc(K * sizeof(*A_u));
			for(int i = 0; i < NUM; ++i)
				time[i] = malloc(K * sizeof(**time));
			first = 0;
		}
		else
		{
			if(sscanf(line, "%lf\t%u", &X_d[c], &A_u[c]) != 2)
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
			naive_int_exp(X_d[i], A_u[i]);
		time[0][i] = clock() - t;

		t = clock();
		for(j=0; j < I; j++)
			squaring_int_exp(X_d[i], A_u[i]);
		time[1][i] = clock() - t;
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
