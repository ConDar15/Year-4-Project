#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>
#include <stdlib.h>

#include "exact_root.h"
#include "newton_sqrt.h"
#include "test_utilities.h"

int main(int argc, char **argv)
{
	FILE *in;
	int first = 1, c = 0, I;
	long unsigned int msec;
	unsigned int K;
	size_t r;
	char line[1024];
	clock_t *time[3], t;

	double *N_d, *T_d;
	//Assumes that a single argument is given and that it is a file of
	//	test values.
	if(argc == 3)
	{
		sscanf(argv[2], "%d", &I);
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
				T_d = malloc(K * sizeof(*T_d));
				time[0] = malloc(K * sizeof(**time));
				time[1] = malloc(K * sizeof(**time));
				time[2] = malloc(K * sizeof(**time));
				first = 0;
			}
			else
			{
				if(sscanf(line, "%lf\t%lf", 
					  &N_d[c], &T_d[c]) != 2)
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
				newton_sqrt_v1(N_d[i], T_d[i]);
			time[0][i] = clock() - t;
			
			t = clock();
			for(j=0; j < I; j++)
				newton_sqrt_v2(N_d[i], T_d[i]);
			time[1][i] = clock() - t;
			
			t = clock();
			for(j=0; j < I; j++)
				newton_sqrt_v3(N_d[i], T_d[i]);
			time[2][i] = clock() - t;
		}
		for(int i = 0; i < 3; i++)
		{
			t = time_average(time[i], K); 
			msec = (long unsigned int)t * 1000 / CLOCKS_PER_SEC;
			printf("Average time for newton_sqrt_v%d: %lu.%04lus\n",
			 	   i+1, msec/1000, msec % 1000);
			t = time_total(time[i], K); 
			msec = (long unsigned int)t * 1000 / CLOCKS_PER_SEC;
			printf("Total time for newton_sqrt_v%d: %lu.%04lus\n",
				   i+1, msec/1000, msec % 1000);
			t = time_min(time[i], K); 
			msec = (long unsigned int)t * 1000 / CLOCKS_PER_SEC;
			printf("Minimum time for newton_sqrt_v%d: %lu.%04lus\n",
				   i+1, msec/1000, msec % 1000);
			t = time_max(time[i], K); 
			msec = (long unsigned int)t * 1000 / CLOCKS_PER_SEC;
			printf("Maximum time for newton_sqrt_v%d: %lu.%04lus\n",
				   i+1, msec/1000, msec % 1000);
			printf("\n");
		}
	}
	else
		printf("Invalid arguments provided\n"
		       "This program expects a single filename as input");

}
