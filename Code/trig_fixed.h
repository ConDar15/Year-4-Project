#ifndef TRIG_FIXED
#define TRIG_FIXED
	#include <inttypes.h>

	#ifdef TEST_FIXED
		#if TEST_FIXED == 8
			#define USE_FIXED_08
		#elif TEST_FIXED == 16
			#define USE_FIXED_16
		#elif TEST_FIXED == 32
			#define USE_FIXED_32
		#elif TEST_FIXED == 64
			#define USE_FIXED_64
		#else
			#error "TEST_FIXED must be defined with a value " \
				   "in {8, 16, 32, 64}\n"
		#endif
	#else
		#if defined __INT64_TYPE__
			#define USE_FIXED_64
		#elif defined __INT32_TYPE__
			#define USE_FIXED_32
		#elif defined __INT16_TYPE__
			#define USE_FIXED_16
		#elif defined __INT8_TYPE__
			#define USE_FIXED_08
		#else
			#error "Available integer types on the system are not " \
				   "defined in system constants\n"
		#endif
	#endif

	#if defined USE_FIXED_08	
		#define TRIG_FIXED_TYPE int8_t
		#define TRIG_FIXED_UTYPE uint8_t
		#define BITS 8
		#define FRAC_BITS 5
		#define UFRAC_BITS 6
		#define CONVERSION_VALUE 0x20
	#elif defined USE_FIXED_16
		#define TRIG_FIXED_TYPE int16_t
		#define TRIG_FIXED_UTYPE uint16_t
		#define BITS 16
		#define FRAC_BIT 13
		#define UFRAC_BITS 14
		#define CONVERSION_VALUE 0x2000
	#elif defined USE_FIXED_32
		#define TRIG_FIXED_TYPE int32_t
		#define TRIG_FIXED_UTYPE uint32_t
		#define BITS 32
		#define FRAC_BIT 29
		#define UFRAC_BITS 30
		#define CONVERSION_VALUE 0x20000000
	#elif defined USE_FIXED_64
		#define TRIG_FIXED_TYPE int64_t
		#define TRIG_FIXED_UTYPE uint64_t
		#define BITS 64
		#define FRAC_BIT 61
		#define UFRAC_BITS 62
		#define CONVERSION_VALUE 0x2000000000000000
	#else
		#error "How did you get here?"
	#endif
#endif
