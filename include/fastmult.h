/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 4/3/09
--------------------------------------------------------------------*/

#ifndef _FASTMULT_H_
#define _FASTMULT_H_

#include <util.h>
#include <dd.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_FHT_POWER 28

#define HUGE_TWIDDLE_CUTOFF 20

typedef struct {
        double *small;
        double *large;
} huge_twiddle_t;

typedef struct {
	uint32 log2_runlength;

	uint32 precision_changed;
	dd_precision_t old_precision;

	volatile double round_constant[2];

	double *twiddle[HUGE_TWIDDLE_CUTOFF + 1];

	huge_twiddle_t huge_twiddle[MAX_FHT_POWER + 1 - HUGE_TWIDDLE_CUTOFF];
} fastmult_info_t;

void fastmult_info_init(fastmult_info_t *info);
void fastmult_info_free(fastmult_info_t *info);

void fastmult(uint32 *a, uint32 awords, 
		uint32 *b, uint32 bwords,
		uint32 *prod, fastmult_info_t *info);

#ifdef __cplusplus
}
#endif

#endif /* _FASTMULT_H_ */
