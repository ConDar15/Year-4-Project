#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "cordic_trig.h"
#include "utilities.h"

//FIXED_ANGLES and FIXED_K_VALUES are defined in "cordic_trig.h",
//	they represent arrays of fixed point pre calculated constants
const cordic_fixed_t TRIG_ANGLES[]  = FIXED_ANGLES;
const cordic_fixed_t TRIG_K_VALUES[] = FIXED_K_VALUES;

double *cordic_trig(const double theta, const unsigned int iter)
{
	//Checks that -pi/2 <= theta <= pi/2 before continuing
	assert(-HALF_PI <= theta && theta <= HALF_PI);
	//If the given iter value is more thant the value of MAX_ITER, or
	//	iter is 0 then n is set to MAX_ITER, otherwise n is set to iter
	unsigned int n = (iter > MAX_ITER || iter == 0) ? MAX_ITER : iter;
	
	//Initialises the initial values for 
	cordic_fixed_t x = TRIG_K_VALUES[n-1], y = 0, t;
	cordic_fixed_t beta = double_to_fixed(theta);
	//Assigns memory for the return array
	double *result = malloc(2*sizeof(*result));

	//Checks for the extreme values of -pi/2, 0 and pi/2
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
	
	//Executes if beta != 0
	if(beta)
	{
		//Main loop
		for(int i = 0; i < n; i++)
		{
			//t is a temporary variable for x, to help in the update
			//	process
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

	result[0] = fixed_to_double(x);
	result[1] = fixed_to_double(y);
	return result;
}

double cordic_atan_bounded(const double z, const unsigned int iter)
{
	//Ensures that the given value is in the range of 0 <= z <= 1
	assert(0 <= z && z <= 1);
	//If the given iter value is more thant the value of MAX_ITER, or
	//	iter is 0 then n is set to MAX_ITER, otherwise n is set to iter
	unsigned int n = (iter > MAX_ITER || iter == 0) ? MAX_ITER : iter;

	//Sets the initial values, note that y/x = z is the important factor
	cordic_fixed_t x = FIXED_ONE >> 1, y = double_to_fixed(z) >> 1, t;
	cordic_fixed_t beta = 0;

	//Checks for the extreme case
	if(y == 0)
		return 0;
	
	//Main loop
	for(int i = 0; i < n; i++)
	{
		t = x;
		if(y < 0)
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
	
	return fixed_to_double(beta);
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

double cordic_acos(double x, unsigned int n)
{
	assert(-1 <= x && x <= 1);
	return x == 0 ? HALF_PI 
				  : x > 0 ? cordic_atan(sqrt(1 - x*x)/x, n)
				          : HALF_PI + cordic_asin(-x, n);
}

double cordic_asin(double x, unsigned int n)
{
	assert(-1 <= x && x <= 1);
	return x == 1 ? HALF_PI
				  : x == -1 ? -HALF_PI
				  			: cordic_atan(x/sqrt(1 - x*x), n);
}

double cordic_atan(double x, unsigned int n)
{
	if(x < 0)
		return -cordic_atan(-x, n);
	
	if(x >= 1)
		return HALF_PI/2 + cordic_atan_bounded((x-1)/(x+1), n);
	
	return cordic_atan_bounded(x, n);
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
			
			case 'd':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("aTan(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_atan(x, n));
				else
					printf("Uasge: %s d <x=value for aTan(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;
			
			case 'e':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("aCos(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_acos(x, n));
				else
					printf("Uasge: %s d <x=value for aCos(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;

			case 'f':
				if(argc == 5 &&
				   sscanf(argv[2], "%lf", &x) == 1 &&
				   sscanf(argv[3], "%u" , &n) == 1 &&
				   sscanf(argv[4], "%u" , &D) == 1)
					printf("aSin(%.*lf) = %.*lf\n",
							d(D), x, D, cordic_asin(x, n));
				else
					printf("Uasge: %s d <x=value for aSin(x)> <n> "
						   "<D=Number of digits to display>\n",
						   argv[0]);
				break;
			
			default:
				printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
		}
	}
	else
		printf("Usage: %s <a/b/c/d/e/f> <arguments>\n", argv[0]);
}
#endif
