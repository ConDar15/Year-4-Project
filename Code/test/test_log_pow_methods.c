#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "taylor_exp_log.h"
#include "test_utilities.h"
#include "hyperbolic_log.h"
#include "cont_frac_exp.h"
#include "int_exp.h"

int main(int argc, char **argv)
{
	FILE *in;
	int first = 1, c = 0, NUM = 6;
	uintmax_t msec;
	unsigned int K, I;
	size_t r;
	char line[1024];
	clock_t *time[6], t;
	char *name[] = {"taylor_log", "hyperbolic_log", "builtin_log", 
					"taylor_pow", "improved_pow", "builtin_pow"};

	double *X_d, *Y_d;
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
			Y_d = malloc(K * sizeof(*Y_d));
			for(int i = 0; i < NUM; ++i)
				time[i] = malloc(K * sizeof(**time));
			first = 0;
		}
		else
		{
			if(sscanf(line, "%lf\t%lf", &X_d[c], &Y_d[c]) != 2)
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
			taylor_log(X_d[i], X_d[i], 35);
		time[0][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			hyperbolic_log(X_d[i], Y_d[i], 8);
		time[1][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			log(Y_d[i])/log(X_d[i]);
		time[2][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			taylor_pow(X_d[i], Y_d[i], 35);
		time[3][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			improved_pow(X_d[i], Y_d[i], 12);
		time[4][i] = clock() - t;
		
		t = clock();
		for(j=0; j < I; j++)
			pow(X_d[i], Y_d[i]);
		time[5][i] = clock() - t;
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
