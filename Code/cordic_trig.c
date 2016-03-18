#include <stdio.h>
#include <stdlib.h>

#include "cordic_trig.h"
#include "utilities.h"

const cordic_fixed_t TRIG_ANGLES[]  = FIXED_ANGLES;
const cordic_fixed_t TRIG_K_VALUES[] = FIXED_K_VALUES;

double *cordic_trig(const double theta, const unsigned int iter)
{
	unsigned int n = (iter > MAX_ITER || iter == 0) ? MAX_ITER : iter;
	
	cordic_fixed_t x = TRIG_K_VALUES[n-1], y = 0, t;
	cordic_fixed_t beta = double_to_fixed(theta);
	
	double *result = malloc(2*sizeof(*result));

	if(beta == 0)
	{
		x = FIXED_ONE;
		beta = 0;
	}
	else if(beta == FIXED_HALF_PI)
	{
		x = 0;
		y = FIXED_ONE;
		beta = 0;
	}
	else if(beta == -FIXED_HALF_PI)
	{
		x = 0;
		y = -FIXED_ONE;
		beta = 0;
	}
	
	if(beta)
	{
		for(int i = 0; i < n; i++)
		{
			t = x;
			if(beta >= 0)
			{
				x = x - (y >> i);
				y = y + (t >> i);
				beta -= TRIG_ANGLES[i];
			}
			else
			{
				x = x + (y >> i);
				y = y - (t >> i);
				beta += TRIG_ANGLES[i];
			}
		}
	}

	result[0] = (double)x / 0x4000000000000000;
	result[1] = (double)y / 0x4000000000000000;
	return result;
}

double cordic_cos(double x, unsigned int n)
{
	double *r;

	if(x >= 0)
	{
		while(x >= TWO_PI)
			x -= TWO_PI;

		if(x >= PI)
			if(x - PI >= HALF_PI)
				r = cordic_trig(TWO_PI - x, n);
			else
			{
				r = cordic_trig(x - PI, n);
				r[0] = -r[0];
			}
		else
			if(x >= HALF_PI)
			{
				r = cordic_trig(PI - x, n);
				r[0] = -r[0];
			}
			else
				r = cordic_trig(x, n);
		return r[0];
	}
	return cordic_cos(-x, n);
}

double cordic_sin(double x, unsigned int n)
{
	return cordic_cos(x - HALF_PI, n);
}

double cordic_tan(double x, unsigned int n)
{
	double *r;
	if(x >= 0)
	{
		while(x >= PI)
			x -= PI;
		
		if(x >= HALF_PI)
		{
			r = cordic_trig(PI - x, n);
			return -1 * r[1]/r[0];
		}
		
		r = cordic_trig(x, n);
		return r[1]/r[0];
	}
	return -1 * cordic_tan(-x, n);
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	double x;
	unsigned int n, D;

	if(argc > 1)
	{
		switch(argv[1][0])
		{
			case 'a':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Cos(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_cos(x, n));
				else
					printf("Uasge: %s a <x=value for Cos(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;
			
			case 'b':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Sin(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_sin(x, n));
				else
					printf("Uasge: %s a <x=value for Sin(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;
			
			case 'c':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("Tan(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_tan(x, n));
				else
					printf("Uasge: %s a <x=value for Tan(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;

			default:
				printf("Usage: %s <a/b/c> <arguments>\n", argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c> <arguments>\n", argv[0]);
}
#endif
