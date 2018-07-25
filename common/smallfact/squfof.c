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

#include <common.h>

/* This SQUFOF implementation will race as many multipliers
   of the input as will fit into 62 bits */

static const uint32 multipliers[] = {1, 3, 5, 7, 
				11, 3*5, 3*7, 3*11, 
				5*7, 5*11, 7*11, 
				3*5*7, 3*5*11, 3*7*11, 
				5*7*11, 3*5*7*11 };

#define MAX_MULTIPLIERS sizeof(multipliers)/sizeof(uint32)
#define QSIZE 50

/* the maximum number of inner loop iterations for all 
   multipliers combined */

#define MAX_CYCLES 100000

/* the number of iterations to do before switching to the
   next multiplier */

#define ONE_CYCLE_ITER 300

typedef struct {
	uint64 n;
	uint32 sqrtn[MAX_MULTIPLIERS]; 
	uint32 cutoff[MAX_MULTIPLIERS];
	uint32 q0[MAX_MULTIPLIERS];
	uint32 p1[MAX_MULTIPLIERS];
	uint32 q1[MAX_MULTIPLIERS];
	uint32 num_saved[MAX_MULTIPLIERS];
	uint16 saved[MAX_MULTIPLIERS][QSIZE];
	uint8 failed[MAX_MULTIPLIERS];
} squfof_data_t;

/*-----------------------------------------------------------------------*/
uint32 squfof_one_cycle(squfof_data_t *data, uint32 mult_idx, 
			uint32 num_iter, uint32 *factor) {

	/* Perform one unit of work for SQUFOF with one multiplier */

	uint32 sqrtn = data->sqrtn[mult_idx];
	uint32 cutoff = data->cutoff[mult_idx];
	uint32 num_saved = data->num_saved[mult_idx];
	uint32 multiplier = 2 * multipliers[mult_idx];
	uint32 coarse_cutoff = cutoff * multiplier;
	uint16 *saved = data->saved[mult_idx];
	uint32 q0 = data->q0[mult_idx];
	uint32 p1 = data->p1[mult_idx];
	uint32 q1 = data->q1[mult_idx];

	uint32 i, j;
	uint32 p0 = 0;
	uint32 sqrtq = 0;
	uint64 scaledn;

	for (i = 0; i < num_iter; i++) {		
		uint32 q, tmp;

		/* compute (sqrtn+p1)/q1; since this is unity
		   more than half the time, special case the
		   division to save some effort */

		tmp = sqrtn + p1 - q1;
		q = 1;
		if (tmp >= q1)
			q += tmp / q1;

		/* compute the next numerator and denominator */
		
		p0 = q * q1 - p1;
		q0 = q0 + (p1 - p0) * q;

		/* In order to avoid trivial factorizations,
		   save q1 values that are small in a list. Only
		   values that do not contain factors of the
		   multiplier must be saved; however, the GCD
		   is 5x more expensive than all the rest of the
		   work performed in the loop. To avoid most GCDs,
		   only attempt one if q1 is less than the cutoff
		   times the full multiplier 
		   
		   If the queue overflows (very rare), signal that
		   this multiplier failed and move on to another one */

		if (q1 < coarse_cutoff) {
			tmp = q1 / mp_gcd_1(q1, multiplier);

			if (tmp < cutoff) {
				if (num_saved >= QSIZE) {
					data->failed[mult_idx] = 1;
					return i;
				}
				saved[num_saved++] = (uint16)tmp;
			}
		}

		/* if q0 is a perfect square, then the factorization
		   has probably succeeded. Most of the squareness
		   tests out there require multiple divisions and 
		   complicated loops. We approximate these tests by
		   ensuring that the last 3 bits are 001 after
		   dividing (shifting) out all factors of 4. */

		tmp = q0;
		while(!(tmp & 3)) {
			tmp >>= 2;
		}
		if ((tmp & 7) == 1) {

			/* q0 is probably a perfect square. Take the
			   square root by cheating */

			sqrtq = (uint32)(sqrt((double)q0));

			if (sqrtq * sqrtq == q0) {

				/* it *is* a perfect square. If it has
				   not appeared previously in the list
				   for this multiplier, then we're almost
				   finished */

				for (j = 0; j < num_saved; j++) {
					if (saved[j] == sqrtq)
						break;
				}

				if (j == num_saved)
					break;
			}
		}

		/* perform the odd half of the SQUFOF cycle */

		tmp = sqrtn + p0 - q0;
		q = 1;
		if (tmp >= q0)
			q += tmp / q0;

		p1 = q * q0 - p0;
		q1 = q1 + (p0 - p1) * q;

		if (q0 < coarse_cutoff) {
			tmp = q0 / mp_gcd_1(q0, multiplier);

			if (tmp < cutoff) {
				if (num_saved >= QSIZE) {
					data->failed[mult_idx] = 1;
					return i;
				}
				saved[num_saved++] = (uint16)tmp;
			}
		}
	}

	if (sqrtq == 1) {

		/* the above found a trivial factor, so this
		   multiplier has failed */

		data->failed[mult_idx] = 1;
		return i;
	}
	else if (i == num_iter) {
	
		/* no square root found; save the parameters
		   and go on to the next multiplier */

		data->q0[mult_idx] = q0;
		data->p1[mult_idx] = p1;
		data->q1[mult_idx] = q1;
		data->num_saved[mult_idx] = num_saved;
		return i;
	}

	/* square root found; the algorithm cannot fail now.
	   Compute the inverse quadratic form and iterate */

	q0 = sqrtq;
	p1 = p0 + sqrtq * ((sqrtn - p0) / sqrtq);
	scaledn = data->n * (uint64)multipliers[mult_idx];
	q1 = (scaledn - (uint64)p1 * (uint64)p1) / (uint64)q0;

	while (1) {
		uint32 q, tmp;

		tmp = sqrtn + p1 - q1;
		q = 1;
		if (tmp >= q1)
			q += tmp / q1;

		p0 = q * q1 - p1;
		q0 = q0 + (p1 - p0) * q;

		if (p0 == p1) {
			q0 = q1;
			break;
		}

		tmp = sqrtn + p0 - q0;
		q = 1;
		if (tmp >= q0)
			q += tmp / q0;

		p1 = q * q0 - p0;
		q1 = q1 + (p0 - p1) * q;

		if (p0 == p1)
			break;
	}

	/* q0 is the factor of n. Remove factors that exist
	   in the multiplier and save whatever's left */

	q0 = q0 / mp_gcd_1(q0, multiplier);
	*factor = q0;
	return i;
}

/*-----------------------------------------------------------------------*/
uint32 squfof(mp_t *n) {

	/* Factor a number up to 62 bits in size using SQUFOF.
	   For n the product of two primes, this routine will
	   succeed with very high probability, although the
	   likelihood of failure goes up as n increases in size.
	   Empirically, 62-bit factorizations fail about 5% of the
	   time; for smaller n the failure rate is nearly zero. 
	   
	   If a factor is found, it is returned. If it is a trivial
	   factor, the return value is one. If SQUFOF failed, the
	   return value is zero */

	uint32 num_mult;
	uint32 factor_found = 0;
	mp_t scaledn, sqrt_scaledn;
	uint32 i, num_iter, num_failed;
	squfof_data_t data;

	/* turn n into a uint64 */

	data.n = (uint64)(n->val[1]) << 32 | (uint64)(n->val[0]);

	/* for each multiplier */

	for (i = 0; i < MAX_MULTIPLIERS; i++) {

		/* use the multiplier if the multiplier times n
		   will fit in 62 bits. Because multipliers are
		   initialized in order of increasing size, when
		   one is too big then all the rest are also too big */

		mp_mul_1(n, multipliers[i], &scaledn);
		if (mp_bits(&scaledn) > 62)
			break;

		mp_isqrt(&scaledn, &sqrt_scaledn);

		/* initialize the rest of the fields for this
		   multiplier */

		data.sqrtn[i] = sqrt_scaledn.val[0];
		data.cutoff[i] = (uint32)(sqrt(2.0 * (double)data.sqrtn[i]));
		data.num_saved[i] = 0;
		data.failed[i] = 0;

		data.q0[i] = 1;
		data.p1[i] = data.sqrtn[i];
		data.q1[i] = scaledn.val[0] - data.p1[i] * data.p1[i];

		/* if n is a perfect square, don't run the algorithm;
		   the factorization has already taken place */

		if (data.q1[i] == 0)
			return data.p1[i];
	}
	if (i == 0)
		return 0;

	/* perform a block of work using each multiplier in
	   turn, until our budget of work for factoring n is
	   exhausted */

	num_mult = i;
	num_iter = 0;
	num_failed = 0;
	while (num_iter < MAX_CYCLES) {
		
		/* for each cycle of multipliers, begin with the
		   multiplier that is largest. These have a higher
		   probability of factoring n quickly */

		for (i = num_mult - 1; (int)i >= 0 ; i--) {
			if (data.failed[i])
				continue;

			/* work on this multiplier for a little while */

			num_iter += squfof_one_cycle(&data, i, ONE_CYCLE_ITER,
							&factor_found);

			/* if all multipliers have failed, then SQUFOF
			   has failed */

			if (data.failed[i]) {
				if (++num_failed == num_mult)
					return 0;
			}
			if (factor_found)
				return factor_found;
		}
	}

	return 0;
}
