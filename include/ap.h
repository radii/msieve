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

#ifndef _AP_H_
#define _AP_H_

#include <mp.h>
#include <fastmult.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Basic arbitrary-precision arithmetic implementation.
   We do not fix the amount of precision, because the
   latency of arithmetic like multiplies and mods is assumed
   to cost far more than the occaisional malloc() or free() */

typedef struct {
	uint32 num_alloc : 31;  /* number of words allocated */
	uint32 sign : 1;	/* POSITIVE or NEGATIVE */
	uint32 nwords;		/* number of nonzero words in val[] */
	uint32 *val;
} ap_t;

	/* initialize an ap_t */

#define ap_init(a) memset(a, 0, sizeof(ap_t))

static INLINE void ap_clear(ap_t *a) {
	if (a->val)
		free(a->val);
	ap_init(a);
}

void ap_copy(ap_t *src, ap_t *dest);

void ap_mp2ap(mp_t *src, uint32 sign, ap_t *dest);

void ap_si2ap(uint32 i, uint32 sign, ap_t *dest);

	/* return the number of bits needed to hold an ap_t.
   	   This is equivalent to floor(log2(a)) + 1. */

uint32 ap_bits(ap_t *a);

	/* Addition and subtraction; a + b = sum
	   or a - b = diff. sum or diff may overwrite 
	   either input operand */

void ap_add(ap_t *a, ap_t *b, ap_t *sum);
void ap_sub(ap_t *a, ap_t *b, ap_t *diff);

	/* return -1, 0, or 1 if a is less than, equal to,
	   or greater than b, respectively (unsigned compare) */

static INLINE int32 ap_cmp_abs(const ap_t *a, const ap_t *b) {

	int32 i;

	if (a->nwords > b->nwords)
		return 1;
	if (a->nwords < b->nwords)
		return -1;

	for (i = a->nwords - 1; i >= 0; i--) {
		if (a->val[i] > b->val[i])
			return 1;
		if (a->val[i] < b->val[i])
			return -1;
	}

	return 0;
}

	/* return -1, 0, or 1 if a is less than, equal to,
	   or greater than b, respectively (signed compare) */

static INLINE int32 ap_cmp(const ap_t *a, const ap_t *b) {

	if (a->sign == POSITIVE && b->sign == NEGATIVE)
		return 1;
	if (a->sign == NEGATIVE && b->sign == POSITIVE)
		return -1;
	
	if (a->sign == POSITIVE)
		return ap_cmp_abs(a, b);
	else
		return -ap_cmp_abs(a, b);
}


	/* quick test for zero or one ap_t */

#define ap_is_zero(a) ((a)->nwords == 0)
#define ap_is_one(a) ((a)->nwords == 1 && (a)->val[0] == 1)

	/* Shift 'a' left by 'shift' bit positions.
	   The result may overwrite 'a' */

void ap_lshift(ap_t *a, uint32 shift, ap_t *res);

	/* Shift 'a' right by 'shift' bit positions.
	   The result may overwrite 'a' */

void ap_rshift(ap_t *a, uint32 shift, ap_t *res);

	/* multiply a by b */

void ap_mul(ap_t *a, ap_t *b, ap_t *prod, fastmult_info_t *info);

	/* compute the integer reciprocal of 'a'. The result
	   is suitable for division operations where the numerator
	   has size <= (div_bits + bits(a)) */

void ap_recip(ap_t *a, ap_t *res, uint32 div_bits, fastmult_info_t *info);

	/* Compute the remainder of num / den. If the 
	   reciprocal value passed in is NULL, a temporary 
	   reciprocal is created. den is assumed positive,
	   and the result is only assured to be less than
	   den in absolute value */

void ap_mod(ap_t *num, ap_t *den, 
		ap_t *recip, ap_t *rem, 
		fastmult_info_t *info);

#define FFT_MIN_WORDS 64

#ifdef __cplusplus
}
#endif

#endif /* !_AP_H_ */
