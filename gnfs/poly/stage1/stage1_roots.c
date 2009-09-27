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

#if 1
#define CHECK
#endif

#if MAX_P_FACTORS > 4
#error "too many factors"
#endif

#define PRIME_P_LIMIT 0xfffff000

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

	free_prime_sieve(&s->prime_sieve);

	free(s->sieve_block);
	free(s->good_primes.primes);
	free(s->bad_primes.primes);
	
	mpz_clear(s->p);
	mpz_clear(s->p2);
	mpz_clear(s->m0);
	mpz_clear(s->nmodp2);
	mpz_clear(s->tmp1);
	mpz_clear(s->tmp2);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_clear(s->accum[i]);
	for (i = 0; i < MAX_ROOTS; i++)
		mpz_clear(s->roots[i]);
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
static uint32
get_prime_roots(poly_search_t *poly, uint32 p, uint32 *roots)
{
	mp_poly_t tmp_poly;
	mp_t *low_coeff;
	uint32 high_coeff;
	uint32 degree = poly->degree;

	memset(&tmp_poly, 0, sizeof(mp_poly_t));
	tmp_poly.degree = degree;
	tmp_poly.coeff[degree].num.nwords = 1;
	tmp_poly.coeff[degree].num.val[0] = p - 1;

	if (mp_gcd_1(p, (uint32)mpz_tdiv_ui(
			poly->high_coeff, (mp_limb_t)p)) > 1)
		return 0;

	low_coeff = &tmp_poly.coeff[0].num;
	low_coeff->val[0] = mpz_tdiv_ui(poly->trans_N, (mp_limb_t)p);
	if (low_coeff->val[0])
		low_coeff->nwords = 1;

	return poly_get_zeros(roots, &tmp_poly, 
				p, &high_coeff, 0);
}

/*------------------------------------------------------------------------*/
void
sieve_fb_init(sieve_fb_t *s, poly_search_t *poly,
		uint32 factor_min, uint32 factor_max)
{
	uint32 i;
	prime_sieve_t prime_sieve;
	uint32 num_squarefree;

	if (factor_max <= factor_min)
		return;

	memset(s, 0, sizeof(sieve_fb_t));

	mpz_init(s->p);
	mpz_init(s->p2);
	mpz_init(s->m0);
	mpz_init(s->nmodp2);
	mpz_init(s->tmp1);
	mpz_init(s->tmp2);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_init(s->accum[i]);
	for (i = 0; i < MAX_ROOTS; i++)
		mpz_init(s->roots[i]);

	s->sieve_block = (uint8 *)xmalloc(SIEVE_SIZE * sizeof(uint8));
	s->degree = poly->degree;
	s->good_primes.num_primes = 0;
	s->good_primes.num_primes_alloc = 500;
	s->good_primes.primes = (sieve_prime_t *)xmalloc(
					s->good_primes.num_primes_alloc * 
					sizeof(sieve_prime_t));
	s->bad_primes.num_primes = 0;
	s->bad_primes.num_primes_alloc = 100;
	s->bad_primes.primes = (sieve_prime_t *)xmalloc(
					s->bad_primes.num_primes_alloc * 
					sizeof(sieve_prime_t));

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
			uint32 roots[MAX_POLY_DEGREE];
			uint32 num_roots = get_prime_roots(poly, p, roots);

			if (num_roots > 0) {
				uint32 logval = (uint32)(LOG_SCALE * 
							log((double)p) / 
							M_LN2 + 0.5);
				sieve_prime_t *sp;

				sieve_add_prime(&s->good_primes, p, logval);

				sp = s->good_primes.primes +
				     s->good_primes.num_primes - 1;
				sp->num_roots = num_roots;
				for (i = 0; i < num_roots; i++)
					sp->roots[i] = roots[i];

			}
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
}

/*------------------------------------------------------------------------*/
static void
sieve_run(sieve_fb_t *s)
{
	uint32 i;
	double cutoff = floor(LOG_SCALE * 
			log((double)s->p_min) / M_LN2 + 0.5);
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
sieve_fb_reset(sieve_fb_t *s, uint64 p_min, uint64 p_max,
		uint32 num_roots_min, uint32 num_roots_max)
{
	uint32 i;
	sieve_prime_list_t *list;

	if (p_min % 2)
		p_min--;
	s->p_min = p_min;
	s->p_max = p_max;
	s->next_prime_p = 0;
	s->next_composite_p = 0;
	s->curr_offset = 0;
	s->num_roots_min = num_roots_min;
	s->num_roots_max = num_roots_max;

	free_prime_sieve(&s->prime_sieve);
	if (p_min >= PRIME_P_LIMIT ||
	    num_roots_min > s->degree) {
		s->allow_prime_p = 0;
	}
	else {
		s->allow_prime_p = 1;
		init_prime_sieve(&s->prime_sieve, (uint32)p_min,
				MIN(p_max, PRIME_P_LIMIT));
	}
	
	list = &s->good_primes;
	for (i = 0; i < list->num_primes; i++) {

		sieve_prime_t *curr = list->primes + i;
		uint32 p = curr->p;
		uint32 rem = p - p_min % p;

		if (rem != p && rem % 2 == 0)
			rem += p;
		curr->r = rem / 2;
	}

	list = &s->bad_primes;
	for (i = 0; i < list->num_primes; i++) {

		sieve_prime_t *curr = list->primes + i;
		uint32 p = curr->p;
		uint32 rem = p - p_min % p;

		if (rem != p && rem % 2 == 0)
			rem += p;
		curr->r = rem / 2;
	}

	sieve_run(s);
}

/*------------------------------------------------------------------------*/
static uint32 
lift_roots(sieve_fb_t *s, poly_search_t *poly, 
		uint64 p, uint32 num_roots)
{
	uint32 i;
	uint32 degree = s->degree;

	uint64_2gmp(p, s->p);
	mpz_mul(s->p2, s->p, s->p);
	mpz_tdiv_r(s->nmodp2, poly->trans_N, s->p2);
	mpz_tdiv_r(s->m0, poly->trans_m0, s->p2);

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
		mpz_add(s->tmp1, poly->trans_m0, s->roots[i]);
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
get_composite_roots_core(sieve_fb_t *s, poly_search_t *poly,
			uint64 p, uint32 num_factors, 
			uint32 *factors)
{
	uint32 i, j, k, i0, i1, i2, i3;
	uint32 crt_p[MAX_P_FACTORS];
	uint32 num_roots[MAX_P_FACTORS];
	uint64 prod[MAX_P_FACTORS];
	uint32 roots[MAX_P_FACTORS][MAX_POLY_DEGREE];
	sieve_prime_t *primes = s->good_primes.primes;
	uint32 degree = s->degree;

	for (i = j = 0; j < MAX_P_FACTORS && i < num_factors; i++, j++) {
		sieve_prime_t *sp = primes + factors[i];
		uint32 power_limit = (uint32)(-1) / sp->p;

		crt_p[j] = sp->p;
		num_roots[j] = sp->num_roots;
		for (k = 0; k < num_roots[j]; k++) {
			roots[j][k] = sp->roots[k];
		}

		while (i < num_factors - 1 && factors[i] == factors[i+1]) {

			uint32 nmodp, new_power;

			if (crt_p[j] > power_limit)
				return 0;

			new_power = crt_p[j] * sp->p;
			nmodp = mpz_tdiv_ui(poly->trans_N, 
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

		return num_roots[0];
	}

	for (i = 0; i < num_factors; i++) {
		prod[i] = p / crt_p[i];
		prod[i] = prod[i] * mp_modinv_1((uint32)(prod[i] %
						crt_p[i]), crt_p[i]);
	}
	mpz_set_ui(s->accum[i], (mp_limb_t)0);
	uint64_2gmp(p, s->p);

	i0 = i1 = i2 = i3 = i = 0;
	switch (num_factors) {
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
		}}}}
	}

	return i;
}

/*------------------------------------------------------------------------*/
#define MAX_FACTORS 20

static uint32
get_composite_roots(sieve_fb_t *s, poly_search_t *poly, uint64 p)
{
	uint32 i, j;
	uint64 psave = p;
	uint32 factors[MAX_FACTORS];
	sieve_prime_t *primes = s->good_primes.primes;
	uint32 num_primes = s->good_primes.num_primes;
	uint32 num_roots;

	for (i = j = 0, num_roots = 1; i < num_primes; i++) {
		uint32 x = primes[i].p;
		if (p % x == 0) {
			num_roots *= primes[i].num_roots;
			do {
				p /= x;
				factors[j++] = i;
			} while (p % x == 0);

			if (p < x)
				break;
		}
	}

	if (p > 1 ||
	    num_roots < s->num_roots_min ||
	    num_roots > s->num_roots_max) {
		return 0;
	}

	return get_composite_roots_core(s, poly, psave, j, factors);
}

/*------------------------------------------------------------------------*/
static uint64
get_next_composite(sieve_fb_t *s)
{
	while (1) {
		uint32 curr_offset;
		uint8 *sieve_block = s->sieve_block;

		for (curr_offset = s->curr_offset; 
				curr_offset < SIEVE_SIZE; curr_offset++) {

			if (!(sieve_block[curr_offset] & 0x80)) 
				continue;

			s->curr_offset = curr_offset+1;
			return s->p_min + (2 * curr_offset + 1);
		}

		sieve_run(s);
		s->p_min += 2 * SIEVE_SIZE;
		s->curr_offset = curr_offset = 0;
		if (s->p_min > s->p_max)
			break;
	}

	return s->p_min;
}

/*------------------------------------------------------------------------*/
uint64
sieve_fb_next(sieve_fb_t *s, poly_search_t *poly,
		root_callback callback, void *extra)
{
	uint32 i;
	uint64 p;
	uint32 num_roots;

	while (1) {
		if (s->allow_prime_p && s->next_prime_p == 0) {
			s->next_prime_p = get_next_prime(&s->prime_sieve);

			if (s->next_prime_p >= PRIME_P_LIMIT)
				s->allow_prime_p = 0;
		}

		if (s->next_composite_p == 0)
			s->next_composite_p = get_next_composite(s);

		if (s->allow_prime_p &&
		    s->next_prime_p < s->next_composite_p) {

			uint32 roots[MAX_POLY_DEGREE];

			p = s->next_prime_p;
			s->next_prime_p = 0;
			num_roots = get_prime_roots(poly, (uint32)p, roots);
			for (i = 0; i < num_roots; i++)
				mpz_set_ui(s->roots[i], (mp_limb_t)roots[i]);
		}
		else {
			p = s->next_composite_p;
			s->next_composite_p = 0;
			if (s->next_prime_p == p)
				s->next_prime_p = 0;

			num_roots = get_composite_roots(s, poly, p);
		}

		if (p > s->p_max)
			break;

		if (num_roots > 0) {
			lift_roots(s, poly, p, num_roots);
			callback(p, num_roots, s->roots, extra);
			break;
		}
	}

	return p;
}

