/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

 $Id$
--------------------------------------------------------------------*/

#ifndef _MP_INT_H_
#define _MP_INT_H_

#include <mp.h>

#ifdef __cplusplus
extern "C" {
#endif
	/* routines that are used within the multiple precision
	   code but are also useful by themselves and as building
	   blocks in other libraries */

	/* return the index of the first nonzero word in
	   x, searching backwards from max_words */

static INLINE uint32 num_nonzero_words(uint32 *x, uint32 max_words) {

	uint32 i;
	for (i = max_words; i && !x[i-1]; i--)
		;
	return i;
}

	/* Internal structure used by routines that need
	   to do division */

typedef struct {
	uint32 nwords;
	uint32 val[2 * MAX_MP_WORDS];
} big_mp_t;

	/* internal division routine; the input is twice the
	   size of an mp_t to allow for a quotient and remainder
	   each that large. Note that no check is made that
	   the quotient fits in an mp_t */

void mp_divrem_core(big_mp_t *num, mp_t *denom, mp_t *quot, mp_t *rem);

#ifdef __cplusplus
}
#endif

#endif /* !_MP_INT_H_ */
