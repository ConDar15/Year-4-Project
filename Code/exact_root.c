#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "utilities.h"
#include "exact_root.h"

char *root_digits_precise(char *N, unsigned int D)
{
	//Counter variables
	unsigned int i, a;
	//The offset value used to set the correct character's value
	unsigned int o = 0;
	//Real and Integer types from GMP and MPFR used for precision
	mpfr_t Yr, Nr,   T, tmpr_0, tmpr_1;
	mpz_t  P,  tmpz, Yz;
	
	//Allocates memory for the number of digits requested plus 5 to be safe
	char *R = malloc((D+5) * sizeof(*R));

	//Sets Nr from the provided string representing N
	mpfr_init_set_str(Nr, N, 10, MPFR_RNDN);
	mpfr_init(Yr);
	mpfr_init(tmpr_0);
	mpfr_init(tmpr_1);
	
	//P will be used to keep track of the current partial solution
	mpz_init_set_ui(P, 0);
	mpz_init(Yz);
	mpz_init(tmpz);

	//T represents the power of the 10 that the current digit represents
	//T is initially floor(n/2) where N == K*10^n and K is in [0, 10)
	//T is of the form 10^t
	mpfr_init(T);
	
	//The mpfr_log10 function is used to help find the exponent (power 10)
	//	of Nr
	mpfr_log10(T, Nr, MPFR_RNDN);
	mpfr_div_ui(T, T, 2, MPFR_RNDN);
	mpfr_floor(T, T);
	//Similar to log, but for exponentiation
	mpfr_exp10(T, T, MPFR_RNDN);

	//This takes into account numbers of the form 0.x
	if(mpfr_cmp_ui(T, 1) < 0)
	{
		R[0] = '0';
		R[1] = '.';
		//Offset set to 2 to indicate there are 2 pre-assigned characters
		o = 2;
	}

	//Main loop
	for(i=0; i <= D; i++)
	{
		//Calculates 10^(2t) and 20*P
		mpfr_mul(tmpr_0, T, T, MPFR_RNDN);
		mpz_mul_ui(tmpz, P, 20);
		//tmpr_1 is used to prevent re-calculation later on
		mpfr_set_ui(tmpr_1, 0, MPFR_RNDN);

		/*
		This loop stops when any digit produces a Y too large or all
			digits have been conisdered.
		In both cases a will be one greater than required and thus
			must be decremented afterwards 
		*/
		for(a=1; a <= 9; a++)
		{
			//Calculates N - (20*P + a)*a*10^(2t)
			mpz_add_ui(Yz, tmpz, a);
			mpz_mul_ui(Yz, Yz, a);
			mpfr_mul_z(Yr, tmpr_0, Yz, MPFR_RNDN);
			
			if(mpfr_cmp(Yr, Nr) > 0)
				//The exit condition for the loop has been met
				break;
			else
				//tmpr_1 updated to remove the need for re-calculation
				mpfr_set(tmpr_1, Yr, MPFR_RNDN);
		}

		//Decrements a and adds the correct digit to the result string
		R[i+o] = DIGITS[--a];
		
		//Reduces Nr by the largest Yr found in the previous loop
		mpfr_sub(Nr, Nr, tmpr_1, MPFR_RNDN);

		//Break if an exact solution is found
		//Note that due to the representation of floating point numbers it
		//	is possible to have found an exact solution with a positive
		//	remainder that is very close to zero. Unfortunately there is
		//	no way to test for this without knowing, the exact precision
		//	of the input beforehand.
		if(mpfr_cmp_ui(Nr, 0) == 0)
		{
			//This loop adds 0s to a string where an exact solution has
			//	been found but needs right padding with zeros.
			while(mpfr_cmp_ui(T, 1) > 0)
			{
				R[++i + o] = '0';
				mpfr_div_ui(T, T, 10, MPFR_RNDN);
			}
			break;
		}
	
		//Calculates P = 10*P + a
		mpz_mul_ui(P, P, 10);
		mpz_add_ui(P, P, a);
		//Calculates T = T/10 => 10^t -> 10^(t-1)
		mpfr_div_ui(T, T, 10, MPFR_RNDN);

		//If we have dropped below 10^0 for the first time then add
		//	a '.' to the result string and increase the offset to 1
		//This case only occurs if no '.' is in the string already
		if(o == 0 && mpfr_cmp_ui(T, 1) < 0)
		{
			o = 1;
			R[i+o] = '.';
		}
	}

	//Adds a null character to terminate the string
	R[i+o+1] = '\0';
	return R;
}

//The use of uintmax_t gives the largest number of unsigned integers 
//	for which this function will work with.
uintmax_t uint_sqrt(uintmax_t num)
{
	//Represents the value of 2r(2^m), where r is the
	//	current known part of the integer root
	uintmax_t res = 0;
	//Represents the largest power of (2^m)^2 = 4^m, the initial value 
	//	is calculated as 011...11 XOR 0011...11 as the size
	//	of uintmax_t is not known beforehand
	uintmax_t bit = (UINTMAX_MAX >> 1) ^ (UINTMAX_MAX >> 2);

	//Finds the largest power of 4 that is at most 'num' in value
	while(bit > num)
		bit >>= 2;

	//while(bit) is equivalent to while(bit > 0) for unsigned integers
	while(bit)
	{
		//Checks the two cases for updating 'res' and 'num'
		if(num >= res + bit)
		{
			//'num' is used to keep track of the difference betweek
			//	r and the original value, N, that was to be rooted.
			num -= res + bit;
			//This calculates 'res' -> 2(r + 2^m)*2^(m-1) using addition 
			//	and	bitshifting 
			res = (res >> 1) + bit;
		}
		//In the other case 'res' -> 2r(2^m-1)
		else
			res >>= 1;
		
		//Move on to the next lower power of 2
		bit >>= 2;
	}

	//Returns the integer part of the square root
	return res;
}

#ifdef COMPILE_MAIN
int main(int argc, char **argv)
{
	uintmax_t N;
	unsigned int p, d;
	char *R;
	
	if (argc == 1)
	{
		printf("Usage: %s [a/b] [arguments]", argv[0]);
		exit(1);
	}
	
	switch(argv[1][0])
	{
		case 'a':
			if(argc == 5 &&
				sscanf(argv[3], "%u", &d) == 1 &&
				sscanf(argv[4], "%u", &p) == 1)
			{
				mpfr_set_default_prec(p);
				printf("sqrt(%s) ~=\n\t%s", argv[2], 
					root_digits_precise(argv[2], d));
			}
			else
				printf("Usage: %s a <N=Number to sqrt> "
					   "<d=Number of significant digits> "
					   "<p=bits of precision to use>\n", argv[0]);
			break;

		case 'b':
			if(argc == 3 &&
				sscanf(argv[2], "%" SCNuMAX, &N) == 1)
				printf("int_sqrt(%" PRIuMAX ") ~= %" PRIuMAX "\n", 
						N, uint_sqrt(N));
			else
				printf("Usage: %s b <N=Positive integer to sqrt>\n",
						argv[0]);
			break;

		default:
			printf("Usage: %s [a/b] [arguments]", argv[0]);
	}
}
#endif
