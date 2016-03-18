#ifndef CORDIC_TRIG_HEADER
#define CORDIC_TRIG_HEADER

#include "trig_fixed.h"
#include "trig_utilities.h"

typedef TRIG_FIXED_TYPE cordic_fixed_t;

#if BITS == 64
	#define FIXED_ONE 0x4000000000000000
	#define NEG_CONSTANT 0x8000000000000000
	#define FIXED_HALF_PI 0x6487ed5110b4611a
	#define FIXED_ANGLES {0x3243f6a8885a308d, 0x1dac670561bb4f68, \
						  0x0fadbafc96406eb1, 0x07f56ea6ab0bdb71, \
					      0x03feab76e59fbd38, 0x01ffd55bba97624a, \
						  0x00fffaaadddb94d5, 0x007fff5556eeea5c, \
						  0x003fffeaaab7776e, 0x001ffffd5555bbbb, \
						  0x000fffffaaaaaddd, 0x0007fffff555556e, \
						  0x0003fffffeaaaaab, 0x0001ffffffd55555, \
						  0x0000fffffffaaaaa, 0x00007fffffff5555, \
						  0x00003fffffffeaaa, 0x00001ffffffffd55, \
						  0x00000fffffffffaa, 0x000007fffffffff5, \
						  0x000003fffffffffe, 0x000001ffffffffff, \
						  0x0000010000000000, 0x0000008000000000, \
						  0x0000004000000000, 0x0000002000000000, \
						  0x0000001000000000, 0x0000000800000000, \
						  0x0000000400000000, 0x0000000200000000, \
						  0x0000000100000000, 0x0000000080000000, \
						  0x0000000040000000, 0x0000000020000000, \
						  0x0000000010000000, 0x0000000008000000, \
						  0x0000000004000000, 0x0000000002000000, \
						  0x0000000001000000, 0x0000000000800000, \
						  0x0000000000400000, 0x0000000000200000, \
						  0x0000000000100000, 0x0000000000080000, \
						  0x0000000000040000, 0x0000000000020000, \
						  0x0000000000010000, 0x0000000000008000, \
						  0x0000000000004000, 0x0000000000002000, \
						  0x0000000000001000, 0x0000000000000800, \
						  0x0000000000000400, 0x0000000000000200, \
						  0x0000000000000100, 0x0000000000000080, \
						  0x0000000000000040, 0x0000000000000020, \
						  0x0000000000000010, 0x0000000000000008, \
						  0x0000000000000004, 0x0000000000000002, \
						  0x0000000000000001}
	#define FIXED_K_VALUES {0x2d413cccfe779921, 0x287a26c490921db6, \
							0x2744c374daf46d2f, 0x26f72283bd67fbda, \
							0x26e3b58305ddeb19, 0x26ded9f57b2c3e7a, \
							0x26dda30d3e4fd185, 0x26dd5552e1641def, \
							0x26dd41e4454da117, 0x26dd3d089dfa47c8, \
							0x26dd3bd1b42095ce, 0x26dd3b83f9a9db95, \
							0x26dd3b708b0c282b, 0x26dd3b6baf64bb03, \
							0x26dd3b6a787adfb4, 0x26dd3b6a2ac068e0, \
							0x26dd3b6a1751cb2b, 0x26dd3b6a127623be, \
							0x26dd3b6a113f39e3, 0x26dd3b6a10f17f6c, \
							0x26dd3b6a10de10ce, 0x26dd3b6a10d93527, \
							0x26dd3b6a10d7fe3d, 0x26dd3b6a10d7b082, \
							0x26dd3b6a10d79d14, 0x26dd3b6a10d79838, \
							0x26dd3b6a10d79701, 0x26dd3b6a10d796b3, \
							0x26dd3b6a10d796a0, 0x26dd3b6a10d7969b, \
							0x26dd3b6a10d7969a, 0x26dd3b6a10d7969a, \
							0x26dd3b6a10d7969a, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699, 0x26dd3b6a10d79699, \
							0x26dd3b6a10d79699}
	#define MAX_ITER 63
#elif BITS == 32
	#define FIXED_ONE 0x40000000
	#define NEG_CONSTANT 0x80000000
	#define FIXED_HALF_PI 0x6487ed51
	#define FIXED_ANGLES {0x3243f6a8, 0x1dac6705, 0x0fadbafc, \
						  0x07f56ea6, 0x03feab76, 0x01ffd55b, \
						  0x00fffaaa, 0x007fff55, 0x003fffea, \
						  0x001ffffd, 0x00100000, 0x00080000, \
						  0x00040000, 0x00020000, 0x00010000, \
						  0x00008000, 0x00004000, 0x00002000, \
						  0x00001000, 0x00000800, 0x00000400, \
						  0x00000200, 0x00000100, 0x00000080, \
						  0x00000040, 0x00000020, 0x00000010, \
						  0x00000008, 0x00000004, 0x00000002, \
						  0x00000001}
	#define FIXED_K_VALUES {0x2d413ccc, 0x287a26c4, 0x2744c374, \
							0x26f72283, 0x26e3b583, 0x26ded9f5, \
							0x26dda30d, 0x26dd5552, 0x26dd41e4, \
							0x26dd3d08, 0x26dd3bd1, 0x26dd3b83, \
							0x26dd3b70, 0x26dd3b6b, 0x26dd3b6a, \
							0x26dd3b6a, 0x26dd3b6a, 0x26dd3b6a, \
							0x26dd3b6a, 0x26dd3b6a, 0x26dd3b6a, \
							0x26dd3b6a, 0x26dd3b6a, 0x26dd3b6a, \
							0x26dd3b6a, 0x26dd3b6a, 0x26dd3b6a, \
							0x26dd3b6a, 0x26dd3b6a, 0x26dd3b6a, \
							0x26dd3b6a}
	#define max_iter 30
#elif BITS == 16
	#define FIXED_ONE 0x4000
	#define NEG_CONSTANT 0x8000
	#define FIXED_HALF_PI 0x6487
	#define FIXED_ANGLES {0x3243, 0x1dac, 0x0fad, 0x07f5, 0x03fe, \
						  0x01ff, 0x0100, 0x0080, 0x0040, 0x0020, \
						  0x0010, 0x0008, 0x0004, 0x0002, 0x0001} 
	#define FIXED_K_VALUES {0x2d41, 0x287a, 0x2744, 0x26f7, 0x26e3, \
							0x26de, 0x26dd, 0x26dd, 0x26dd, 0x26dd, \
							0x26dd, 0x26dd, 0x26dd, 0x26dd, 0x26dd}
	#define MAX_ITER 15
#elif BITS == 8
	#define FIXED_ONE 0x40
	#define NEG_CONSTANT 0x80
	#define FIXED_HALF_PI 0x64
	#define FIXED_ANGLES {0x32, 0x1d, 0x10, 0x08, 0x04, 0x02, 0x01}
	#define FIXED_K_VALUES {0x2d, 0x28, 0x27, 0x26, 0x26, 0x26, 0x26}	
	#define MAX_ITER 7
#else
	#error "you shouldn't be able to get here; you done messed up"
#endif

	double *cordic_trig(const double, const unsigned int);
	
	double cordic_cos(double, unsigned int);
	double cordic_sin(double, unsigned int);
	double cordic_tan(double, unsigned int);
	
	double cordic_atan_bounde(const double, const unsigned int);
	double cordic_atan(double, unsigned int);
	double cordic_acos(double, unsigned int);
	double cordic_asin(double, unsigned int);

#endif
