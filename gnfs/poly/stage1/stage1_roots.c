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

#include "stage1.h"

#define SIEVE_SIZE 16384
#define LOG_SCALE 3.5

#if 0
#define CHECK
#endif

#if MAX_P_FACTORS > 6
#error "too many factors"
#endif

/*------------------------------------------------------------------------*/
static uint32
lift_root_32(uint32 n, uint32 r, uint32 old_power, 
		uint32 p, uint32 d)
{
	uint32 q;
	uint32 p2 = old_power * p;
	uint64 rsave = r;

	q = mp_modsub_1(n % p2, mp_expo_1(r, d, p2), p2) / old_power;
	r = mp_modmul_1(d, mp_expo_1(r % p, d - 1, p), p);
	r = mp_modmul_1(q, mp_modinv_1(r, p), p);
	return rsave + old_power * r;
}

/*------------------------------------------------------------------------*/
void
sieve_fb_free(sieve_fb_t *s)
{
	uint32 i;

	free(s->sieve_block);
	free(s->good_primes.primes);
	free(s->bad_primes.primes);
	free(s->small_roots);

	for (i = 0; i < s->max_roots; i++)
		mpz_init(s->roots[i]);
	free(s->roots);

	mpz_clear(s->p);
	mpz_clear(s->p2);
	mpz_clear(s->nmodp2);
	mpz_clear(s->m0);
	mpz_clear(s->tmp1);
	mpz_clear(s->tmp2);
	for (i = 0; i < MAX_P_FACTORS + 1; i++)
		mpz_clear(s->accum[i]);
}

/*------------------------------------------------------------------------*/
static void
sieve_add_prime(sieve_prime_list_t *list, uint32 p, uint32 logval)
{
	sieve_prime_t *curr;

	if (list->num_primes == list->num_primes_alloc) {
		list->num_primes_alloc *= 2;
		list->primes = (sieve_prime_t *)xrealloc(list->primes,
						list->num_primes_alloc * 
						sizeof(sieve_prime_t));	
	}
	curr = list->primes + list->num_primes++;
	curr->p = p;
	curr->log_p = (uint8)logval;
}

/*------------------------------------------------------------------------*/
void
sieve_fb_init(sieve_fb_t *s, poly_batch_t *poly,
		uint32 factor_min, uint32 factor_max)
{
	uint32 i, j, k;
	mp_poly_t tmp_poly;
	prime_sieve_t prime_sieve;
	uint32 num_squarefree;
	uint32 max_roots;
	uint32 degree = poly->degree;

	memset(s, 0, sizeof(sieve_fb_t));

	if (factor_max <= factor_min)
		return;

	for (i = 1, max_roots = degree; i < MAX_P_FACTORS; i++)
		max_roots *= degree;

	s->max_roots = max_roots;
	s->roots = (mpz_t *)xmalloc(max_roots * sizeof(mpz_t));
	for (i = 0; i < max_roots; i++)
		mpz_init(s->roots[i]);
	mpz_init(s->p);
	mpz_init(s->p2);
	mpz_init(s->nmodp2);
	mpz_init(s->m0);
	mpz_init(s->tmp1);
	mpz_init(s->tmp2);
	for (i = 0; i < MAX_P_FACTORS + 1; i++)
		mpz_init(s->accum[i]);

	s->sieve_block = (uint8 *)xmalloc(SIEVE_SIZE * sizeof(uint8));
	s->good_primes.num_primes = 0;
	s->good_primes.num_primes_alloc = 500;
	s->good_primes.primes = (sieve_prime_t *)xmalloc(
					s->good_primes.num_primes_alloc * 
					sizeof(sieve_prime_t));
	s->bad_primes.num_primes = 0;
	s->bad_primes.num_primes_alloc = 10;
	s->bad_primes.primes = (sieve_prime_t *)xmalloc(
					s->bad_primes.num_primes_alloc * 
					sizeof(sieve_prime_t));
	s->num_small_roots = 0;
	s->num_small_roots_alloc = 1000;
	s->small_roots = (uint32 *)xmalloc(s->num_small_roots_alloc * 
						sizeof(uint32));

	init_prime_sieve(&prime_sieve, 2, factor_max);

	while (1) {
		uint32 p = get_next_prime(&prime_sieve);

		if (p >= factor_max) {
			break;
		}
		else if (p <= factor_min) {
			sieve_add_prime(&s->bad_primes, p, 0);
			continue;
		}
		else {
			uint32 logval = (uint32)(LOG_SCALE * log((double)p) / 
							M_LN2 + 0.5);
			sieve_add_prime(&s->good_primes, p, logval);
		}
	}

	free_prime_sieve(&prime_sieve);
	num_squarefree = s->good_primes.num_primes;

	for (i = 0; i < num_squarefree; i++) {

		uint32 p = s->good_primes.primes[i].p;
		uint32 power_limit = factor_max / p;

		if (power_limit < p) {
			break;
		}
		else {
			uint32 power = p;
			uint32 logval = (uint32)(LOG_SCALE * log((double)p) / 
							M_LN2 + 0.5);
			while (power < power_limit) {
				power *= p;
				sieve_add_prime(&s->good_primes, power, logval);
			}
		}
	}

	memset(&tmp_poly, 0, sizeof(mp_poly_t));
	tmp_poly.degree = degree;
	tmp_poly.coeff[degree].num.nwords = 1;

	for (i = 0; i < num_squarefree; i++) {

		sieve_prime_t *curr_p = s->good_primes.primes + i;
		uint32 p = curr_p->p;

		curr_p->max_roots = 0;
		tmp_poly.coeff[degree].num.val[0] = p - 1;

		for (j = 0; j < poly->num_poly; j++) {

			uint32 roots[MAX_POLY_DEGREE];
			uint32 high_coeff;
			uint32 num_roots = 0; 
			curr_poly_t *curr_poly = poly->batch + j;

			if (mp_gcd_1(p, (uint32)mpz_tdiv_ui(
					curr_poly->high_coeff, 
					(mp_limb_t)p)) == 1) {

				mp_t *low_coeff = &tmp_poly.coeff[0].num;
				low_coeff->val[0] = mpz_tdiv_ui(
							curr_poly->trans_N, 
							(mp_limb_t)p);
				if (low_coeff->val[0])
					low_coeff->nwords = 1;

				num_roots = poly_get_zeros(roots, &tmp_poly, 
							p, &high_coeff, 0);
			}

			curr_p->root_offsets[j] = s->num_small_roots;
			if (num_roots == 0)
				continue;

			curr_p->max_roots = MAX(curr_p->max_roots, num_roots);

			if (s->num_small_roots + num_roots >= 
					s->num_small_roots_alloc) {
				s->num_small_roots_alloc = 2 * 
						(s->num_small_roots + 
						 num_roots);
				s->small_roots = (uint32 *)xrealloc(
						s->small_roots,
						s->num_small_roots_alloc *
						sizeof(uint32));
			}
			for (k = 0; k < num_roots; k++)
				s->small_roots[s->num_small_roots+k] = roots[k];
			s->num_small_roots += k;
		}

		curr_p->root_offsets[j] = s->num_small_roots;
	}
}

/*------------------------------------------------------------------------*/
static void
sieve_run(sieve_fb_t *s)
{
	uint32 i;
	double cutoff = floor(LOG_SCALE * 
			log((double)s->base) / M_LN2 + 0.5);
	uint8 *sieve_block = s->sieve_block;
	uint32 num_good_primes = s->good_primes.num_primes;
	uint32 num_bad_primes = s->bad_primes.num_primes;
	sieve_prime_t *good_primes = s->good_primes.primes;
	sieve_prime_t *bad_primes = s->bad_primes.primes;

	memset(sieve_block, (int)(cutoff - 2 * LOG_SCALE), 
			(size_t)SIEVE_SIZE);

	for (i = 0; i < num_good_primes; i++) {
		sieve_prime_t *curr = good_primes + i;
		uint32 p = curr->p;
		uint32 r = curr->r;
		uint8 log_p = curr->log_p;

		while (r < SIEVE_SIZE) {
			sieve_block[r] -= log_p;
			r += p;
		}
		curr->r = r - SIEVE_SIZE;
	}

	for (i = 0; i < num_bad_primes; i++) {
		sieve_prime_t *curr = bad_primes + i;
		uint32 p = curr->p;
		uint32 r = curr->r;

		while (r < SIEVE_SIZE) {
			sieve_block[r] = 0;
			r += p;
		}
		curr->r = r - SIEVE_SIZE;
	}

}

/*------------------------------------------------------------------------*/
void 
sieve_fb_reset(sieve_fb_t *s, uint64 base)
{
	uint32 i;
	sieve_prime_list_t *list;

	if (base % 2)
		base--;
	s->base = base;
	s->curr_offset = 0;

	list = &s->good_primes;
	for (i = 0; i < list->num_primes; i++) {

		sieve_prime_t *curr = list->primes + i;
		uint32 p = curr->p;
		uint32 rem = p - base % p;

		if (rem != p && rem % 2 == 0)
			rem += p;
		curr->r = rem / 2;
	}

	list = &s->bad_primes;
	for (i = 0; i < list->num_primes; i++) {

		sieve_prime_t *curr = list->primes + i;
		uint32 p = curr->p;
		uint32 rem = p - base % p;

		if (rem != p && rem % 2 == 0)
			rem += p;
		curr->r = rem / 2;
	}

	sieve_run(s);
}

/*------------------------------------------------------------------------*/
static uint32 
lift_roots(sieve_fb_t *s, curr_poly_t *curr, 
		uint64 p, uint32 num_roots,
		uint32 degree)
{
	uint32 i;

	uint64_2gmp(p, s->p);
	mpz_mul(s->p2, s->p, s->p);
	mpz_tdiv_r(s->nmodp2, curr->trans_N, s->p2);
	mpz_tdiv_r(s->m0, curr->trans_m0, s->p2);

	for (i = 0; i < num_roots; i++) {

		mpz_powm_ui(s->tmp1, s->roots[i], (mp_limb_t)degree, s->p2);
		mpz_sub(s->tmp1, s->nmodp2, s->tmp1);
		if (mpz_cmp_ui(s->tmp1, (mp_limb_t)0) < 0)
			mpz_add(s->tmp1, s->tmp1, s->p2);
		mpz_tdiv_q(s->tmp1, s->tmp1, s->p);

		mpz_powm_ui(s->tmp2, s->roots[i], (mp_limb_t)(degree-1), s->p);
		mpz_mul_ui(s->tmp2, s->tmp2, (mp_limb_t)degree);
		mpz_invert(s->tmp2, s->tmp2, s->p);

		mpz_mul(s->tmp1, s->tmp1, s->tmp2);
		mpz_tdiv_r(s->tmp1, s->tmp1, s->p);
		mpz_addmul(s->roots[i], s->tmp1, s->p);
		mpz_sub(s->roots[i], s->roots[i], s->m0);
		if (mpz_cmp_ui(s->roots[i], (mp_limb_t)0) < 0)
			mpz_add(s->roots[i], s->roots[i], s->p2);
#ifdef CHECK
		mpz_add(s->tmp1, curr->trans_m0, s->roots[i]);
		mpz_tdiv_r(s->tmp1, s->tmp1, s->p2);
		mpz_powm_ui(s->tmp1, s->tmp1, (mp_limb_t)degree, s->p2);
		if (mpz_cmp(s->tmp1, s->nmodp2) != 0) {
			gmp_printf("error: %Zd n %Zd p %Zd p2 %Zd m0 %Zd r %Zd\n",
					s->tmp1, s->nmodp2, s->p, s->p2, 
					s->m0, s->roots[i]);
			exit(-1);
		}
#endif
	}

	return num_roots;
}

/*------------------------------------------------------------------------*/
static uint32 
find_composite_roots(sieve_fb_t *s, poly_batch_t *poly,
		uint32 which_poly, uint64 p, 
		uint32 num_factors, uint32 *factors)
{
	uint32 i, j, k, i0, i1, i2, i3, i4, i5;
	uint32 crt_p[MAX_P_FACTORS];
	uint32 num_roots[MAX_P_FACTORS];
	uint64 prod[MAX_P_FACTORS];
	uint32 roots[MAX_P_FACTORS][MAX_POLY_DEGREE];
	curr_poly_t *curr = poly->batch + which_poly;
	sieve_prime_t *primes = s->good_primes.primes;
	uint32 degree = poly->degree;

	for (i = 0; i < num_factors; i++) {
		sieve_prime_t *sp = primes + factors[i];
		uint32 *offset = sp->root_offsets + which_poly;

		if (offset[0] == offset[1])
			return 0;
	}

	for (i = j = 0; j < MAX_P_FACTORS && i < num_factors; i++, j++) {
		sieve_prime_t *sp = primes + factors[i];
		uint32 power_limit = (uint32)(-1) / sp->p;
		uint32 *offset = sp->root_offsets + which_poly;

		crt_p[j] = sp->p;
		num_roots[j] = offset[1] - offset[0];
		for (k = 0; k < num_roots[j]; k++) {
			roots[j][k] = s->small_roots[offset[0] + k];
		}

		while (i < num_factors - 1 && factors[i] == factors[i+1]) {

			uint32 nmodp, new_power;

			if (crt_p[j] > power_limit)
				return 0;

			new_power = crt_p[j] * sp->p;
			nmodp = mpz_tdiv_ui(curr->trans_N, 
						(mp_limb_t)new_power);

			for (k = 0; k < num_roots[j]; k++) {
				roots[j][k] = lift_root_32(nmodp, roots[j][k],
							crt_p[j], sp->p, 
							degree);
			}
			crt_p[j] = new_power;
			i++;
		}
	}
	if (i < num_factors)
		return 0;
	num_factors = j;

	if (num_factors == 1) {
		for (i = 0; i < num_roots[0]; i++)
			mpz_set_ui(s->roots[i], (mp_limb_t)roots[0][i]);

		return lift_roots(s, curr, p, num_roots[0], degree);
	}

	for (i = 0; i < num_factors; i++) {
		prod[i] = p / crt_p[i];
		prod[i] = prod[i] * mp_modinv_1((uint32)(prod[i] %
						crt_p[i]), crt_p[i]);
	}
	mpz_set_ui(s->accum[i], (mp_limb_t)0);
	uint64_2gmp(p, s->p);

	i0 = i1 = i2 = i3 = i4 = i5 = i = 0;
	switch (num_factors) {
	case 6:
		for (i5 = num_roots[5] - 1; (int32)i5 >= 0; i5--) {
			uint64_2gmp(prod[5], s->accum[5]);
			mpz_mul_ui(s->accum[5], s->accum[5], 
						(mp_limb_t)roots[5][i5]);
	case 5:
		for (i4 = num_roots[4] - 1; (int32)i4 >= 0; i4--) {
			uint64_2gmp(prod[4], s->accum[4]);
			mpz_mul_ui(s->accum[4], s->accum[4], 
						(mp_limb_t)roots[4][i4]);
			mpz_add(s->accum[4], s->accum[4], s->accum[5]);
	case 4:
		for (i3 = num_roots[3] - 1; (int32)i3 >= 0; i3--) {
			uint64_2gmp(prod[3], s->accum[3]);
			mpz_mul_ui(s->accum[3], s->accum[3], 
						(mp_limb_t)roots[3][i3]);
			mpz_add(s->accum[3], s->accum[3], s->accum[4]);
	case 3:
		for (i2 = num_roots[2] - 1; (int32)i2 >= 0; i2--) {
			uint64_2gmp(prod[2], s->accum[2]);
			mpz_mul_ui(s->accum[2], s->accum[2], 
						(mp_limb_t)roots[2][i2]);
			mpz_add(s->accum[2], s->accum[2], s->accum[3]);
	case 2:
		for (i1 = num_roots[1] - 1; (int32)i1 >= 0; i1--) {
			uint64_2gmp(prod[1], s->accum[1]);
			mpz_mul_ui(s->accum[1], s->accum[1], 
						(mp_limb_t)roots[1][i1]);
			mpz_add(s->accum[1], s->accum[1], s->accum[2]);

		for (i0 = num_roots[0] - 1; (int32)i0 >= 0; i0--) {
			uint64_2gmp(prod[0], s->accum[0]);
			mpz_mul_ui(s->accum[0], s->accum[0], 
						(mp_limb_t)roots[0][i0]);
			mpz_add(s->accum[0], s->accum[0], s->accum[1]);

			mpz_tdiv_r(s->accum[0], s->accum[0], s->p);
			mpz_set(s->roots[i++], s->accum[0]);
		}}}}}}
	}

	return lift_roots(s, curr, p, i, degree);
}

/*------------------------------------------------------------------------*/
#define MAX_FACTORS 20

uint64
sieve_fb_next(sieve_fb_t *s, poly_batch_t *poly,
		p_batch_t *pbatch, uint64 limit,
		uint32 num_roots_min)
{
	uint32 i, j;
	uint64 p, psave;
	uint32 factors[MAX_FACTORS];
	uint8 *sieve_block = s->sieve_block;
	sieve_prime_t *primes = s->good_primes.primes;
	uint32 num_primes = s->good_primes.num_primes;
	uint32 curr_offset;
	uint32 done;

	while (1) {

		for (curr_offset = s->curr_offset; 
				curr_offset < SIEVE_SIZE; curr_offset++) {

			if (!(sieve_block[curr_offset] & 0x80)) 
				continue;

			psave = p = s->base + (2 * curr_offset + 1);

			if (p > limit) {
				s->curr_offset = curr_offset;
				return p;
			}

			for (i = j = 0; i < num_primes; i++) {
				uint32 x = primes[i].p;
				if (p % x == 0) {
					do {
						p /= x;
						factors[j++] = i;
					} while (p % x == 0);

					if (p < x)
						break;
				}
			}
			if (p > 1)
				continue;

			if (num_roots_min > 1) {
				uint32 num_roots = 1;
				for (i = 0; i < j; i++) {
					sieve_prime_t *sp = primes + factors[i];
					num_roots *= sp->max_roots;
				}
				if (num_roots < num_roots_min)
					continue;
			}

			for (i = done = 0; i < poly->num_poly; i++) {

				uint32 num_lifted = find_composite_roots(s, 
						poly, i, psave, j, factors);

				if (num_lifted) {
					p_batch_add(pbatch, psave, i,
							num_lifted, s->roots);
					done = 1;
				}
			}
			if (done) {
				s->curr_offset = curr_offset+1;
				return psave;
			}
		}

		sieve_run(s);
		s->base += 2 * SIEVE_SIZE;
		s->curr_offset = curr_offset = 0;
		if (s->base > limit)
			return s->base;
	}

	return s->base;
}
