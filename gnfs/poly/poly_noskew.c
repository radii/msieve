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

/* Polynomial generator for the number field sieve.
   The following generates non-skewed polynomials, and does
   not attempt to optimize them.

   Note that just covering the search space of candidate
   polynomials efficiently is a very interesting problem. The
   ideal solution
   	- covers as much of the space as possible
	- as fast as possible
	- without performing redundant work
	- and covers a broad range of polynomials early on
	- while allowing faster CPUs to explore more of the space
	- in a way that parallelizes easily

   The current implementation does all of these, though it is
   not set up to run in parallel. Hopefully some of the techniques
   used can be applied to skewed polynomials */

#include "poly.h"

/* We specialize to 5th degree polynomials for now */
#define POLY_DEGREE 5

/* The search process works by picking the leading coefficient
   of the algebraic poly to be a smooth number that is a fraction
   between MAX_SHRINKAGE and MIN_SHRINKAGE of the common root 
   'm' of the rational and algebraic polynomials. 
   Once a leading algebraic coefficient is chosen, we analyze 
   any polynomial that has all of its remaining coefficients 
   smaller than (1-WORST_SHRINKAGE)*m in size. The 0.02 figure is
   straight out of Murphy's thesis */

#define MAX_SHRINKAGE 0.00001
#define MIN_SHRINKAGE 0.004
#define WORST_SHRINKAGE 0.02

/* leading coefficients are handled in batches 
   for efficiency */
#define COFACTOR_BATCH_SIZE 100

/* selection-algorithm-specific data */

typedef struct {
	dd_t binomial_coeff[6];  /* see below */
	mp_t *cofactor_batch;    /* see below */
} poly_noskew_t;

/* Leading coefficients of the algebraic polynomial are assumed
   to be the product of two groups of factors. The two groups 
   are coprime to each other; one group is a product of powers
   of small primes up to 31 and is found by sieving. The other
   group is enumerated explicitly */

#define COFACTOR_SIEVE_SIZE 4096

/* primes in the first group */

static const uint8 cofactor_p[] = 
	{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

/* we sieve with powers of the primes in cofactor_p. Each 
   power gets the same log value, and log values stack up 
   at locations divisible by powers of a prime */

static const uint8 cofactor_sieve_p[] = 
	{8, 8, 8, 9, 9, 9, 5, 5, 5, 5, 7, 7, 7, 11, 11,
	13, 13, 17, 17, 19, 19, 23, 23, 29, 29, 31, 31};

/* the sieving stride used in each prime power */

static const uint16 cofactor_sieve_powers[] = 
	{8, 64, 512, 9, 81, 729, 5, 25, 125, 625, 
	 7, 49, 343, 11, 121, 13, 169, 17, 289, 
	 19, 361, 23, 529, 29, 841, 31, 961};

#define NUM_COFACTOR_PRIMES sizeof(cofactor_p)
#define NUM_COFACTOR_ENTRIES sizeof(cofactor_sieve_p)


/* main structure controlling cofactor sieving */

typedef struct {
	uint32 strides[NUM_COFACTOR_ENTRIES];	/* for factors in group 1 */
	uint32 roots[NUM_COFACTOR_ENTRIES];
	uint8 logs[NUM_COFACTOR_ENTRIES];
	uint8 cutoff;
	uint8 *sieve;
	
	uint32 num_large_factors;	/* for factors in group 2 */
	uint32 *large_factors;
} cofactor_sieve_t;

static void cofactor_sieve_init(mp_t *base, cofactor_sieve_t *c);
static void cofactor_sieve_free(cofactor_sieve_t *c);
static void cofactor_sieve_update_logs(mp_t *new_base, cofactor_sieve_t *c);
static void cofactor_sieve_next(cofactor_sieve_t *c);
static void init_cofactor_large(cofactor_sieve_t *c, 
				uint32 min_leftover,
				uint32 max_leftover);


static void search_coeff_batch(mp_t *n, 
			poly_config_t *config, 
			poly_noskew_t *noskew_data,
			uint32 *large_factors,
			uint32 min_large_factor,
			uint32 max_large_factor,
			uint32 default_cofactor);

static void check_poly(mp_t *n, poly_config_t *config,
		       mp_t *m, mp_t *high_coeff);

/*------------------------------------------------------------------*/
#define PACKED_SIEVE_MASK ((uint64)0x80808080 << 32 | 0x80808080)
#define MIN_BIAS 100

void find_poly_noskew(msieve_obj *obj, mp_t *n,
			poly_config_t *config,
			uint32 deadline) {

	/* Search for good NFS polynomials among all those
	   whose leading coefficient is within a given range */

	uint32 i, j, k;
	double log_n, min_coeff, max_coeff;
	mp_t mp_min, mp_max;
	uint32 default_cofactor;
	uint32 min_leftover, max_leftover;
	uint32 min_large_factor, max_large_factor;
	cofactor_sieve_t cofactor_sieve;
	double d, scaled_min_coeff, scaled_max_coeff;
	uint32 blocks_done;
	time_t start_time, curr_time;
	poly_noskew_t noskew_data;

	dd_precision_t precision = 0;
	uint32 precision_changed = 0;

	if (!dd_precision_is_ieee()) {
		precision_changed = 1;
		precision = dd_set_precision_ieee();
	}

	/* initialize polynomial-method-specific data */

	noskew_data.binomial_coeff[0] = dd_div_dxd(-1., 5.);
	noskew_data.binomial_coeff[1] = dd_div_dxd(1.*6, 5.*5*2);
	noskew_data.binomial_coeff[2] = dd_div_dxd(-1.*6*11, 5.*5*5*6);
	noskew_data.binomial_coeff[3] = dd_div_dxd(1.*6*11*16, 5.*5*5*5*24);
	noskew_data.binomial_coeff[4] = dd_div_dxd(-1.*6*11*16*21, 
						5.*5*5*5*5*120);
	noskew_data.binomial_coeff[5] = dd_div_dxd(1.*6*11*16*21*26, 
						5.*5*5*5*5*5*720);

	noskew_data.cofactor_batch = (mp_t *)xmalloc(COFACTOR_BATCH_SIZE * 
						sizeof(mp_t));

	/* determine the range of leading polynomial 
	   coefficients that will be searched */

	log_n = mp_log(n);
	min_coeff = exp((POLY_DEGREE * log(MAX_SHRINKAGE) + log_n) / 
						(POLY_DEGREE + 1));
	max_coeff = exp((POLY_DEGREE * log(MIN_SHRINKAGE) + log_n) / 
						(POLY_DEGREE + 1));

	/* all leading coefficients will be divisible by
	   default_cofactor. This is chosen both to ease
	   the sieving burden and to give all polynomials
	   more roots for the really important small primes */

	if (min_coeff < 1e13)
		default_cofactor = 2*2*3;
	else if (min_coeff < 1e15)
		default_cofactor = 2*2*2*3*5;
	else if (min_coeff < 5e17)
		default_cofactor = 2*2*2*3*5*7;
	else if (min_coeff < 1e20)
		default_cofactor = 2*2*2*3*3*5*5*7;
	else
		default_cofactor = 2*2*2*3*3*5*5*7*11;

	scaled_min_coeff = min_coeff / ((double)default_cofactor * MIN_BIAS);
	scaled_max_coeff = max_coeff / default_cofactor;
	min_leftover = cofactor_p[NUM_COFACTOR_PRIMES-1] + 1;
	max_leftover = (uint32)(scaled_max_coeff / scaled_min_coeff);
	mp_d2mp(&scaled_min_coeff, &mp_min);
	mp_d2mp(&scaled_max_coeff, &mp_max);

	/* leading coefficients are the product of a smooth cofactor
	   and a large factor, which are coprime to each other. The large
	   factor starts off between MIN_BIAS and the largest factor
	   that will create a leading coefficient less than max_coeff */

	cofactor_sieve_init(&mp_min, &cofactor_sieve);
	init_cofactor_large(&cofactor_sieve, min_leftover, max_leftover);

	max_large_factor = cofactor_sieve.num_large_factors;
	min_large_factor = max_large_factor;
	for (i = 0; i < max_large_factor; i++) {
		if (cofactor_sieve.large_factors[i] > MIN_BIAS) {
			min_large_factor = i;
			break;
		}
	}

	/* while time has not run out and the search has
	   not been interrupted and there are still valid
	   large factors */

	blocks_done = 0;
	time(&start_time);
	while (min_large_factor < max_large_factor) {

		uint64 *packed_sieve = (uint64 *)cofactor_sieve.sieve;

		/* find a sequence of smooth cofactors by 
		   sieving. Buffer the next COFACTOR_BATCH_SIZE 
		   smooth cofactors. These are very rare,
		   so check multiple sieve values in parallel.

		   An easy way to parallelize this process is
		   to compress the sieve by a factor k if there
		   are k machines searching for polynomials. 
		   Each search machine starts the sieve at a
		   different offset mod k. Note that this is 
		   *not* the same as multiplying the value of
		   default_cofactor by k! */

		i = 0;
		while (i < COFACTOR_BATCH_SIZE) {
			cofactor_sieve_next(&cofactor_sieve);

			for (j = 0; j < COFACTOR_SIEVE_SIZE/8; j += 4) {
				if (((packed_sieve[j+0] | 
				      packed_sieve[j+1] |
				      packed_sieve[j+2] | 
				      packed_sieve[j+3]) & 
				      PACKED_SIEVE_MASK) == (uint64)0)
					continue;

				for (k = 0; k < 32; k++) {
					if (cofactor_sieve.sieve[8*j+k]&0x80) {
						mp_add_1(&mp_min, 8*j+k, 
						  noskew_data.cofactor_batch+i);
						break;
					}
				}
				if (++i == COFACTOR_BATCH_SIZE)
					break;
			}
			mp_add_1(&mp_min, COFACTOR_SIEVE_SIZE, &mp_min);
		}

		/* search the entire batch at the same time */

		search_coeff_batch(n, config, &noskew_data,
				cofactor_sieve.large_factors, 
				min_large_factor, max_large_factor,
				default_cofactor);

		/* since the values sieved have gotten larger, the
		   range of large factors that still produce a leading
		   coefficient between min_coeff and max_coeff must
		   shrink by a little */

		d = mp_mp2d(&mp_min);
		max_leftover = (uint32)(scaled_max_coeff / d);
		while (max_large_factor &&
			cofactor_sieve.large_factors[max_large_factor-1] > 
							max_leftover) {
			max_large_factor--;
		}

		min_leftover = (uint32)(min_coeff / (d * default_cofactor));
		while (min_large_factor &&
			cofactor_sieve.large_factors[min_large_factor-1] > 
							min_leftover) {
			min_large_factor--;
		}

		cofactor_sieve_update_logs(&mp_min, &cofactor_sieve);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

		/* print progress message */

		if (++blocks_done % 20 == 0) {
			time(&curr_time);
			fprintf(stderr, "%u%% done (processed %u blocks)\r", 
				MIN(100, (uint32)(100.0 * difftime(
				    curr_time, start_time) / deadline + 0.5)),
				blocks_done);
			fflush(stderr);
			if (difftime(curr_time, start_time) > deadline)
				break;
		}
	}

	fprintf(stderr, "\n");
	free(noskew_data.cofactor_batch);
	cofactor_sieve_free(&cofactor_sieve);

	if (precision_changed)
		dd_clear_precision(precision);
}

/*------------------------------------------------------------------*/
/* magic constants for rounding floating
   point values to integers. These require IEEE
   53-bit double precision arithmetic. Do *not*
   make this array static; we need to fool the
   compiler into not optimizing away arithmetic
   that appears redundant */

double dround[2] = {
	6755399441055744.0,
	6755399441055744.0,
};

static void search_coeff_batch(mp_t *n, 
			poly_config_t *config, 
			poly_noskew_t *noskew_data,
			uint32 *large_factors,
			uint32 min_large_factor,
			uint32 max_large_factor,
			uint32 default_cofactor) {
	
	/* do the core work of searching through a batch of
	   polynomials. The list of leading coefficients to look
	   through is the cross product of two sets: the batch
	   of smooth cofactors in config->cofactor_batch and 
	   the numbers in large_factors from index min_large_factor
	   to index max_large_factor. Since the two sets will each
	   never have any repeated elements, and members of the
	   two sets are coprime, we are guaranteed that every
	   product of smooth cofactor and large factor has not been
	   previously searched.

	   In addition, every leading coefficient is divisible
	   by default_cofactor */

	uint32 i, j;
	mp_t tmp1, tmp2;
	double round0 = dround[0];
	double round1 = dround[1];
	mp_t *cofactor_batch = noskew_data->cofactor_batch;
	dd_t *binomial_coeff = noskew_data->binomial_coeff;
	mp_t m, coeff[6];
	dd_t prev_m, curr_m; 

	/* for each large factor */

	for (i = min_large_factor; i < max_large_factor; i++) {

		/* make leading coefficients out of all
		   the smooth cofactors in the batch */

		for (j = 0; j < COFACTOR_BATCH_SIZE; j++) {
			double q, recip_dm, max_dm, dm_half, x;
			double dm, d0, d1, d2, d3, d4, d5;
			double c0, c1, c2, c3, c4, c5;
			mp_t high_coeff;
			int32 correction;
	
			/* form the leading coefficient */

			mp_copy(cofactor_batch + j, &high_coeff);
			mp_mul_1(&high_coeff, default_cofactor, &high_coeff);
			mp_mul_1(&high_coeff, large_factors[i], &high_coeff);
	
			/* the value of 'm' to start with is
			   floor((n / leading_coefficient) ^ (1/5)) */

			if (j == 0) {

				/* For the first cofactor in the list, compute
				   this directly */

				mp_div(n, &high_coeff, &tmp1);
				mp_iroot(&tmp1, 5, &m);
				curr_m = prev_m = dd_mp2dd(&m);
			}
			else {
				/* For the rest of the batch, use the fact
				   that the previous value of m is known.
				   For default cofactor D, large factor L and
				   smooth cofactor S, the value of 'm' for
				   polynomial j is (n / (D*S[j]*L))^(1/5). The
				   previous value of 'm' has all numbers the
				   same except it uses S[j-1]. Thus the current
				   value of m is

				   (n / (D*(S[j-1]+x)*L))^(1/5)
				   = (n / (D*L*S[j-1]*(1 + x/S[j-1])))^(1/5)
				   = (previous m) * (1 + x/S[j-1])^(-1/5)

				   The sequence S[] is the result of a sieving
				   procedure, so x above is usually several
				   digits smaller than S[j-1]. The second term
				   above can be expanded as a binomial series 
				   that converges in only a few terms, and 
				   summing the series (in extended-precision 
				   floating point) is much faster than taking 
				   a multiple-precision root, even if the 
				   previous value of m is used as an
				   approximation to start the Newton iteration.
				   The other parts of the search process run 
				   so fast that root extraction would take a 
				   major fraction of the total time.

				   Wow, the first calculus I've used in about
				   a decade */

				dd_t ddq;

				c0 = cofactor_batch[j].val[1] * MP_RADIX +
					cofactor_batch[j].val[0];
				d0 = cofactor_batch[j-1].val[1] * MP_RADIX +
					cofactor_batch[j-1].val[0];
	
				ddq = dd_div_dxd(c0 - d0, d0);
				curr_m = dd_mul_dd(ddq, binomial_coeff[4]);
				curr_m = dd_add_dd(curr_m, binomial_coeff[3]);
				curr_m = dd_mul_dd(ddq, curr_m);
				curr_m = dd_add_dd(curr_m, binomial_coeff[2]);
				curr_m = dd_mul_dd(ddq, curr_m);
				curr_m = dd_add_dd(curr_m, binomial_coeff[1]);
				curr_m = dd_mul_dd(ddq, curr_m);
				curr_m = dd_add_dd(curr_m, binomial_coeff[0]);
				curr_m = dd_mul_dd(ddq, curr_m);
				curr_m = dd_add_d(curr_m, 1.0);
				curr_m = prev_m = dd_mul_dd(curr_m, prev_m);
				dd_dd2mp(curr_m, &m);
			}
	
			/* Now that the mess above computed m, express
			   n in base m and convert each coefficient to 
			   a double. We don't need all coefficients to
			   full precision here, since all we will care 
			   about is the top few digits of each (more 
			   specifically in how big each is) */

			mp_divrem(    n, &m, &tmp1, coeff + 0);
			mp_divrem(&tmp1, &m, &tmp2, coeff + 1);
			mp_divrem(&tmp2, &m, &tmp1, coeff + 2);
			mp_divrem(&tmp1, &m, &tmp2, coeff + 3);
			mp_divrem(&tmp2, &m, coeff + 5, coeff + 4);
	
			dm = curr_m.hi;
			d0 = mp_mp2d(coeff + 0);
			d1 = mp_mp2d(coeff + 1);
			d2 = mp_mp2d(coeff + 2);
			d3 = mp_mp2d(coeff + 3);
			d4 = mp_mp2d(coeff + 4);
			d5 = mp_mp2d(&high_coeff);
	
			/* Here we use Murphy's trick for non-skewed
			   polynomials: given m, we can easily enumerate
			   all the values of x where converting n from
			   base m to base (m+x) would keep the second-
			   highest polynomial coefficient from exceeding max_dm
			   in size. This lets us fix two coefficients out
			   of 6 and only check the other 4 if this preliminary
			   test passes */

			max_dm = WORST_SHRINKAGE * dm;
			dm_half = dm / 2;
			recip_dm = 1.0 / dm;
	
			/* because the computation of m was approximate,
			   m may be off by one. If it is, correct the
			   second-highest coefficient */

			correction = high_coeff.val[0] - coeff[5].val[0];
			d4 -= correction * dm;
			x = -floor((max_dm - d4) / (5 * d5));
	
			/* now convert from base m to base (m+x), reducing
			   each coefficient to be less than m/2 in absolute
			   value. Note that we really should be reducing mod
			   (m+x) and not mod m, but all that matters at this
			   stage are the high-order digits of each coefficient,
			   so we fudge the division. Base conversion stops
			   when a coefficient is encountered that is too
			   large */

			for (c5 = d5; ; x++) {
				c4 = d4 - 5 * x * d5;
				if (c4 > max_dm)
					continue;
				if (c4 < -max_dm)
					break;
		
				c3 = d3 + x * (-4 * d4 + x * 10 * d5);
				q = c3 * recip_dm + round0 - round1;
				c3 -= q * dm;
				if (c3 < -dm_half) c3 += dm;
				if (c3 > dm_half) c3 -= dm;
				if (c3 < -max_dm || c3 > max_dm)
					continue;
		
				c2 = d2 + x * (-3 * d3 + x * 
						(6 * d4 - x * 10 * d5));
				q = c2 * recip_dm + round0 - round1;
				c2 -= q * dm;
				if (c2 < -dm_half) c2 += dm;
				if (c2 > dm_half) c2 -= dm;
				if (c2 < -max_dm || c2 > max_dm)
					continue;
		
				c1 = d1 + x * (-2 * d2 + x * 
					(3 * d3 + x * (-4 * d4 + x * 5 * d5)));
				q = c1 * recip_dm + round0 - round1;
				c1 -= q * dm;
				if (c1 < -dm_half) c1 += dm;
				if (c1 > dm_half) c1 -= dm;
				if (c1 < -max_dm || c1 > max_dm)
					continue;
		
				c0 = d0 + x * (-d1 + x * (d2 + x * 
						(-d3 + x * (d4 - x * d5))));
				q = c0 * recip_dm + round0 - round1;
				c0 -= q * dm;
				if (c0 < -dm_half) c0 += dm;
				if (c0 > dm_half) c0 -= dm;
				if (c0 < -max_dm || c0 > max_dm)
					continue;
		
				if (x < 0)
					mp_sub_1(&m, (uint32)fabs(x), &tmp1);
				else
					mp_add_1(&m, (uint32)fabs(x), &tmp1);
	
				check_poly(n, config, &tmp1, &high_coeff);
			}
		}
	}
}

/*------------------------------------------------------------------*/
static void check_poly(mp_t *n, poly_config_t *config,
		      mp_t *m, mp_t *high_coeff) {

	/* analyze a polynomial for sieving goodness */

	uint32 i;
	mp_t tmp1, tmp2;
	poly_select_t poly;
	mp_poly_t *rpoly;
	mp_poly_t *apoly;
	
	memset(&poly, 0, sizeof(poly));
	rpoly = &poly.rpoly;
	apoly = &poly.apoly;
	poly.skewness = 1.0;

	/* Write out n in base m */

	mp_divrem(    n, m, &tmp1, &apoly->coeff[0].num);
	mp_divrem(&tmp1, m, &tmp2, &apoly->coeff[1].num);
	mp_divrem(&tmp2, m, &tmp1, &apoly->coeff[2].num);
	mp_divrem(&tmp1, m, &tmp2, &apoly->coeff[3].num);
	mp_divrem(&tmp2, m, &apoly->coeff[5].num, &apoly->coeff[4].num);

	/* change the coefficients to be less than m/2
	   in absolute value */

	mp_rshift(m, 1, &tmp1);
	apoly->degree = 5;
	for (i = 0; i < 5; i++) {
		apoly->coeff[i].sign = POSITIVE;
		if (mp_cmp(&apoly->coeff[i].num, &tmp1) > 0) {
			apoly->coeff[i].sign = NEGATIVE;
			mp_sub(m, &apoly->coeff[i].num, 
					&apoly->coeff[i].num);
			mp_add_1(&apoly->coeff[i+1].num, 1, 
					&apoly->coeff[i+1].num);
		}
	}
	apoly->coeff[i].sign = POSITIVE;

	/* check that the leading coefficient is as expected */

	if (mp_cmp(high_coeff, &apoly->coeff[5].num) != 0) {
		printf("warning: skipping corrupt polynomial\n");
		return;
	}

	/* form the rational polynomial */

	rpoly->degree = 1;
	mp_copy(m, &rpoly->coeff[0].num);
	rpoly->coeff[0].sign = NEGATIVE;
	rpoly->coeff[1].num.nwords = 1;
	rpoly->coeff[1].num.val[0] = 1;
	rpoly->coeff[1].sign = POSITIVE;

	/* perform the analysis */

	analyze_poly(config, &poly);
	save_poly(config, &poly);
}

/*------------------------------------------------------------------*/
#define COFACTOR_LOG_SCALE 2.2

static void cofactor_sieve_init(mp_t *base, cofactor_sieve_t *c) {

	/* set up the cofactor sieve. Sieving will 
	   begin from 'base' and proceed in blocks of 
	   size COFACTOR_SIEVE_SIZE */

	uint32 i;
	double dbase = mp_mp2d(base);

	c->sieve = (uint8 *)aligned_malloc((size_t)COFACTOR_SIEVE_SIZE, 16);
	c->cutoff = (uint8)(COFACTOR_LOG_SCALE * log(dbase) / M_LN2 + 0.5);
	
	for (i = 0; i < NUM_COFACTOR_ENTRIES; i++) {
		uint32 root;
		c->strides[i] = cofactor_sieve_powers[i];
		c->logs[i] = (uint8)(COFACTOR_LOG_SCALE * 
				log((double)cofactor_sieve_p[i]) / 
				M_LN2 + 0.5);
		root = mp_mod_1(base, c->strides[i]);
		if (root > 0)
			root = c->strides[i] - root;
		c->roots[i] = root;
	}
	c->cutoff -= c->logs[i-1];
}

static void cofactor_sieve_free(cofactor_sieve_t *c) {
	aligned_free(c->sieve);
	free(c->large_factors);
}

static void cofactor_sieve_update_logs(mp_t *new_base, cofactor_sieve_t *c) {

	/* as sieve values grow, it becomes harder and harder
	   to find a smooth cofactor that is a fixed fraction
	   of their size. Update the sieving cutoff to reflect 
	   this */

	double dbase = mp_mp2d(new_base);

	c->cutoff = (uint8)(COFACTOR_LOG_SCALE * log(dbase) / M_LN2 + 0.5);
	c->cutoff -= c->logs[NUM_COFACTOR_ENTRIES-1];
}

static void cofactor_sieve_next(cofactor_sieve_t *c) {
	uint32 i;
	uint8 *sieve = c->sieve;

	/* fill the next sieve block. Sieve values begin
	   with the target log value and are decremented
	   as primes update the block. Values that meet
	   the cutoff get their sign bit set, which is
	   easy to check in parallel */

	memset(sieve, c->cutoff, (size_t)COFACTOR_SIEVE_SIZE);
	for (i = 0; i < NUM_COFACTOR_ENTRIES; i++) {
		uint32 r = c->roots[i];
		uint32 s = c->strides[i];
		uint8 logp = c->logs[i];

		while (r < COFACTOR_SIEVE_SIZE) {
			sieve[r] -= logp;
			r += s;
		}
		c->roots[i] = r - COFACTOR_SIEVE_SIZE;
	}
}

static void init_cofactor_large(cofactor_sieve_t *c, 
				uint32 min_leftover,
				uint32 max_leftover) {

	/* find all the numbers between min_leftover and
	   max_leftover that do not have any of the primes
	   in cofactor_p as factors */

	uint32 i, j;
	uint32 sieve_size;
	uint8 *sieve;

	if (min_leftover & 1)
		min_leftover--;
	
	sieve_size = max_leftover - min_leftover + 1;
	sieve = (uint8 *)xmalloc((size_t)(sieve_size / 8 + 1));
	memset(sieve, 0, (size_t)(sieve_size / 8 + 1));

	/* do some quick sieving */

	for (i = 1; i < sizeof(cofactor_p); i++) {
		uint32 p = cofactor_p[i];
		uint32 r = p - (min_leftover % p);

		if (r == p)
			r = 0;

		while (r < sieve_size) {
			sieve[r/8] |= 1 << (r % 8);
			r += p;
		}
	}

	/* count and then read in the survivors */

	for (i = 1, j = 0; i < sieve_size; i += 2) {
		if (!(sieve[i/8] & (1 << (i % 8))))
			j++;
	}

	c->num_large_factors = j;
	c->large_factors = (uint32 *)xmalloc(j * sizeof(uint32));
	for (i = 1, j = 0; i < sieve_size; i += 2) {
		if (!(sieve[i/8] & (1 << (i % 8))))
			c->large_factors[j++] = min_leftover + i;
	}
	free(sieve);
}
