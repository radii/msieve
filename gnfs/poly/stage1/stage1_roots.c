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

#include <stage1.h>

#define SIEVE_SIZE 16384
#define LOG_SCALE 2.2

#define P_PRIME_LIMIT 0xfffff000

#define INVALID_NUM_ROOTS ((uint32)(-1))

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

	free_prime_sieve(&s->p_prime);

	free(s->aprog_data.aprogs);

	free(s->p_sieve.sieve_block);
	free(s->p_sieve.good_primes.primes);
	free(s->p_sieve.bad_primes.primes);
	
	free(s->p_enum.curr_list.list);
	free(s->p_enum.new_list.list);

	mpz_clear(s->p);
	mpz_clear(s->p2);
	mpz_clear(s->m0);
	mpz_clear(s->nmodp2);
	mpz_clear(s->tmp1);
	mpz_clear(s->tmp2);
	for (i = 0; i <= MAX_P_FACTORS; i++)
		mpz_clear(s->accum[i]);
	for (i = 0; i < MAX_ROOTS; i++)
		mpz_clear(s->roots[i]);
}

/*------------------------------------------------------------------------*/
static void
sieve_add_prime(sieve_prime_list_t *list, uint32 p)
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
	curr->log_p = (uint8)(LOG_SCALE * log((double)p) / M_LN2 + 0.5);
}

/*------------------------------------------------------------------------*/
static uint32
get_prime_roots(poly_search_t *poly, uint32 which_poly,
		uint32 p, uint32 *roots)
{
	mp_poly_t tmp_poly;
	mp_t *low_coeff;
	uint32 high_coeff;
	uint32 degree = poly->degree;
	curr_poly_t *c = poly->batch + which_poly;

	memset(&tmp_poly, 0, sizeof(mp_poly_t));
	tmp_poly.degree = degree;
	tmp_poly.coeff[degree].num.nwords = 1;
	tmp_poly.coeff[degree].num.val[0] = p - 1;

	if (mp_gcd_1(p, (uint32)mpz_tdiv_ui(
			c->high_coeff, (mp_limb_t)p)) > 1)
		return 0;

	low_coeff = &tmp_poly.coeff[0].num;
	low_coeff->val[0] = mpz_tdiv_ui(c->trans_N, (mp_limb_t)p);
	if (low_coeff->val[0])
		low_coeff->nwords = 1;

	return poly_get_zeros(roots, &tmp_poly, 
				p, &high_coeff, 0);
}

/*------------------------------------------------------------------------*/
static uint32
sieve_add_aprog(sieve_fb_t *s, poly_search_t *poly, 
		uint32 p, uint32 fb_roots_min, 
		uint32 fb_roots_max)
{
	uint32 i, j;
	uint32 found;
	aprog_list_t *list = &s->aprog_data;
	aprog_t *a;

	if (list->num_aprogs == list->num_aprogs_alloc) {
		list->num_aprogs_alloc *= 2;
		list->aprogs = (aprog_t *)xrealloc(list->aprogs,
						list->num_aprogs_alloc *
						sizeof(aprog_t));
	}

	a = list->aprogs + list->num_aprogs;
	a->p = p;
	for (i = found = 0; i < poly->num_poly; i++) {
		uint32 roots[MAX_POLYSELECT_DEGREE];
		uint32 num_roots = get_prime_roots(poly, i, p, roots);

		a->num_roots[i] = num_roots;
		if (num_roots == 0)
			continue;

		if (num_roots < fb_roots_min ||
		    num_roots > fb_roots_max) {
			found = 0;
			break;
		}

		for (j = 0; j < num_roots; j++)
			a->roots[i][j] = roots[j];
		found++;
	}

	if (found)
		list->num_aprogs++;
	return found;
}

/*------------------------------------------------------------------------*/
void
sieve_fb_init(sieve_fb_t *s, poly_search_t *poly,
		uint32 factor_min, uint32 factor_max,
		uint32 fb_roots_min, uint32 fb_roots_max)
{
	uint32 i, j;
	prime_sieve_t prime_sieve;
	aprog_list_t *aprog = &s->aprog_data;
	p_sieve_t *p_sieve = &s->p_sieve;
	p_enum_t *p_enum = &s->p_enum;

	if (factor_max <= factor_min)
		return;

	memset(s, 0, sizeof(sieve_fb_t));
	s->degree = poly->degree;

	mpz_init(s->p);
	mpz_init(s->p2);
	mpz_init(s->m0);
	mpz_init(s->nmodp2);
	mpz_init(s->tmp1);
	mpz_init(s->tmp2);
	for (i = 0; i <= MAX_P_FACTORS; i++)
		mpz_init(s->accum[i]);
	for (i = 0; i < MAX_ROOTS; i++)
		mpz_init(s->roots[i]);

	i = 500;
	aprog->num_aprogs = 0;
	aprog->num_aprogs_alloc = i;
	aprog->aprogs = (aprog_t *)xmalloc(i * sizeof(aprog_t));

	p_sieve->sieve_block = (uint8 *)xmalloc(SIEVE_SIZE * sizeof(uint8));
	p_sieve->good_primes.num_primes = 0;
	p_sieve->good_primes.num_primes_alloc = i;
	p_sieve->good_primes.primes = (sieve_prime_t *)xmalloc(
					i * sizeof(sieve_prime_t));
	p_sieve->bad_primes.num_primes = 0;
	p_sieve->bad_primes.num_primes_alloc = i;
	p_sieve->bad_primes.primes = (sieve_prime_t *)xmalloc(
					i * sizeof(sieve_prime_t));

	p_enum->curr_list.num_entries_alloc = i;
	p_enum->curr_list.list = (ss_t *)xmalloc(i * sizeof(ss_t));
	p_enum->new_list.num_entries_alloc = i;
	p_enum->new_list.list = (ss_t *)xmalloc(i * sizeof(ss_t));

	init_prime_sieve(&prime_sieve, 2, factor_max);

	while (1) {
		uint32 p = get_next_prime(&prime_sieve);

		if (p >= factor_max)
			break;
		else if (p <= factor_min)
			sieve_add_prime(&p_sieve->bad_primes, p);
		else if (!sieve_add_aprog(s, poly, p,
					fb_roots_min, fb_roots_max))
			sieve_add_prime(&p_sieve->bad_primes, p);
		else
			sieve_add_prime(&p_sieve->good_primes, p);
	}

	free_prime_sieve(&prime_sieve);
	j = p_sieve->good_primes.num_primes;
	for (i = 0; i < j; i++) {
		uint32 p = p_sieve->good_primes.primes[i].p;
		uint32 power_limit = factor_max / p;
		uint32 power = p;

		if (power_limit < p)
			break;

		while (power < power_limit) {
			power *= p;
			sieve_add_prime(&p_sieve->good_primes, power);
		}
	}
}

/*------------------------------------------------------------------------*/
static void
sieve_run(sieve_fb_t *s)
{
	uint32 i;
	p_sieve_t *p_sieve = &s->p_sieve;
	double cutoff = floor(LOG_SCALE * 
			log((double)s->p_min) / M_LN2 + 0.5);
	uint8 *sieve_block = p_sieve->sieve_block;
	uint32 num_good_primes = p_sieve->good_primes.num_primes;
	uint32 num_bad_primes = p_sieve->bad_primes.num_primes;
	sieve_prime_t *good_primes = p_sieve->good_primes.primes;
	sieve_prime_t *bad_primes = p_sieve->bad_primes.primes;

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
static void add_to_enum(subset_sum_t *new_list, ss_t new_entry)
{
	if (new_list->num_entries == new_list->num_entries_alloc) {
		new_list->num_entries_alloc *= 2;
		new_list->list = (ss_t *)xrealloc(new_list->list,
					new_list->num_entries_alloc *
					sizeof(ss_t));
	}
	new_list->list[new_list->num_entries++] = new_entry;
}

#define MAX_POWER 3

static void enum_run(sieve_fb_t *s)
{
	uint32 i, j;
	p_enum_t *p_enum = &s->p_enum;
	uint32 curr_entries = p_enum->curr_list.num_entries;
	subset_sum_t *curr_list = &p_enum->curr_list;
	subset_sum_t *new_list = &p_enum->new_list;
	subset_sum_t tmp;

	aprog_t *aprog = s->aprog_data.aprogs + p_enum->next_prime;
	uint32 p = aprog->p;
	float log_p = log((double)p);

	uint32 power_limit = (uint32)(-1) / p;
	uint32 num_powers;

	i = p;
	for (num_powers = 1; num_powers <= MAX_POWER; num_powers++) {
		if (i >= power_limit)
			break;
		i *= p;
	}
			
	new_list->num_entries = 0;
	for (i = 0; i < curr_entries; i++) {

		ss_t new_ss = curr_list->list[i];

		if (new_ss.log_prod > p_enum->log_p_min ||
		    new_ss.num_factors == MAX_P_FACTORS ||
		    new_ss.log_prod + log_p > p_enum->log_p_max)
			continue;

		add_to_enum(new_list, new_ss);

		new_ss.num_factors++;
		for (j = 0; j < num_powers; j++) {
			if (new_ss.log_prod + log_p > p_enum->log_p_max)
				break;

			new_ss.prod *= p;
			new_ss.log_prod += log_p;
			add_to_enum(new_list, new_ss);
		}
	}

	p_enum->next_prime++;
	tmp = p_enum->curr_list;
	p_enum->curr_list = p_enum->new_list;
	p_enum->new_list = tmp;
}

/*------------------------------------------------------------------------*/
void 
sieve_fb_reset(sieve_fb_t *s, uint64 p_min, uint64 p_max,
		uint32 num_roots_min, uint32 num_roots_max)
{
	uint32 i;

	free_prime_sieve(&s->p_prime);

	if (p_min % 2)
		p_min--;
	s->p_min = p_min;
	s->p_max = p_max;
	s->num_roots_min = num_roots_min;
	s->num_roots_max = num_roots_max;

	if (p_max > 10000000) {
		p_enum_t *p_enum = &s->p_enum;

		s->curr_algo = ALGO_ENUM;

		p_enum->log_p_min = log((double)p_min);
		p_enum->log_p_max = log((double)p_max);
		p_enum->next_prime = 0;
		p_enum->curr_entry = 0;

		p_enum->curr_list.num_entries = 1;
		p_enum->curr_list.list[0].num_factors = 0;
		p_enum->curr_list.list[0].log_prod = 0;
		p_enum->curr_list.list[0].prod = 1;
	}
	else {
		p_sieve_t *p_sieve = &s->p_sieve;
		sieve_prime_list_t *list;

		s->curr_algo = ALGO_SIEVE;
		p_sieve->sieve_start = p_min;

		list = &p_sieve->good_primes;
		for (i = 0; i < list->num_primes; i++) {

			sieve_prime_t *curr = list->primes + i;
			uint32 p = curr->p;
			uint32 rem = p - p_min % p;

			if (rem != p && rem % 2 == 0)
				rem += p;
			curr->r = rem / 2;
		}

		list = &p_sieve->bad_primes;
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
}

/*------------------------------------------------------------------------*/
static uint32 
lift_roots(sieve_fb_t *s, curr_poly_t *c, 
		uint64 p, uint32 num_roots)
{
	uint32 i;
	uint32 degree = s->degree;

	uint64_2gmp(p, s->p);
	mpz_mul(s->p2, s->p, s->p);
	mpz_tdiv_r(s->nmodp2, c->trans_N, s->p2);
	mpz_sub(s->tmp1, c->trans_m0, c->mp_sieve_size);
	mpz_tdiv_r(s->m0, s->tmp1, s->p2);

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
	}

	return num_roots;
}

/*------------------------------------------------------------------------*/
static uint32 
get_composite_roots(sieve_fb_t *s, curr_poly_t *c,
			uint32 which_poly, uint64 p, 
			uint32 num_factors, 
			uint32 *factors,
			uint32 num_roots_min,
			uint32 num_roots_max)
{
	uint32 i, j, k, i0, i1, i2, i3, i4, i5, i6;
	uint32 crt_p[MAX_P_FACTORS];
	uint32 num_roots[MAX_P_FACTORS];
	uint64 prod[MAX_P_FACTORS];
	uint32 roots[MAX_P_FACTORS][MAX_POLYSELECT_DEGREE];
	aprog_t *aprogs = s->aprog_data.aprogs;
	uint32 degree = s->degree;

	for (i = 0, j = 1; i < num_factors; i++) {
		aprog_t *a;

		if (i > 0 && factors[i] == factors[i-1])
			continue;

		a = aprogs + factors[i];
		if (a->num_roots[which_poly] == 0)
			return 0;

		j *= a->num_roots[which_poly];
	}
	if (j < num_roots_min || j > num_roots_max)
		return INVALID_NUM_ROOTS;

	for (i = j = 0; j < MAX_P_FACTORS && i < num_factors; i++, j++) {
		aprog_t *a = aprogs + factors[i];
		uint32 power_limit;

		num_roots[j] = a->num_roots[which_poly];
		crt_p[j] = a->p;
		power_limit = (uint32)(-1) / a->p;
		for (k = 0; k < num_roots[j]; k++) {
			roots[j][k] = a->roots[which_poly][k];
		}

		while (i < num_factors - 1 && factors[i] == factors[i+1]) {

			uint32 nmodp, new_power;

			if (crt_p[j] > power_limit)
				return 0;

			new_power = crt_p[j] * a->p;
			nmodp = mpz_tdiv_ui(c->trans_N, (mp_limb_t)new_power);

			for (k = 0; k < num_roots[j]; k++) {
				roots[j][k] = lift_root_32(nmodp, roots[j][k],
							crt_p[j], a->p, 
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

	i0 = i1 = i2 = i3 = i4 = i5 = i6 = i = 0;
	switch (num_factors) {
	case 7:
		for (i6 = num_roots[6] - 1; (int32)i6 >= 0; i6--) {
			uint64_2gmp(prod[6], s->accum[6]);
			mpz_mul_ui(s->accum[6], s->accum[6], 
						(mp_limb_t)roots[6][i6]);
			mpz_add(s->accum[6], s->accum[6], s->accum[7]);
	case 6:
		for (i5 = num_roots[5] - 1; (int32)i5 >= 0; i5--) {
			uint64_2gmp(prod[5], s->accum[5]);
			mpz_mul_ui(s->accum[5], s->accum[5], 
						(mp_limb_t)roots[5][i5]);
			mpz_add(s->accum[5], s->accum[5], s->accum[6]);
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
		}}}}}}}
	}

	return i;
}

/*------------------------------------------------------------------------*/
#define MAX_FACTORS 20

static uint32
get_composite_factors(sieve_fb_t *s, uint64 p,
			uint32 factors[MAX_FACTORS])
{
	uint32 i, j;
	aprog_t *aprogs = s->aprog_data.aprogs;
	uint32 num_primes = s->aprog_data.num_aprogs;

	for (i = j = 0; i < num_primes; i++) {
		uint32 x = aprogs[i].p;
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
		return 0;

	return j;
}

/*------------------------------------------------------------------------*/
static uint64
get_next_composite(sieve_fb_t *s)
{
	while (1) {
		uint32 curr_offset;
		p_sieve_t *p_sieve = &s->p_sieve;
		uint8 *sieve_block = p_sieve->sieve_block;

		for (curr_offset = p_sieve->curr_offset; 
				curr_offset < SIEVE_SIZE; curr_offset++) {

			if (!(sieve_block[curr_offset] & 0x80)) 
				continue;

			p_sieve->curr_offset = curr_offset + 1;
			return p_sieve->sieve_start + (2 * curr_offset + 1);
		}

		sieve_run(s);
		p_sieve->sieve_start += 2 * SIEVE_SIZE;
		p_sieve->curr_offset = curr_offset = 0;
		if (p_sieve->sieve_start > s->p_max)
			break;
	}

	return P_SEARCH_DONE;
}

/*------------------------------------------------------------------------*/
static uint64
get_next_enum_composite(sieve_fb_t *s)
{
	p_enum_t *p_enum = &s->p_enum;

	while (1) {
		uint32 curr_entry;
		ss_t *list = p_enum->curr_list.list;
		uint32 num_entries = p_enum->curr_list.num_entries;

		for (curr_entry = p_enum->curr_entry;
				curr_entry < num_entries; curr_entry++) {

			ss_t *curr_prod = list + curr_entry;

			if (curr_prod->log_prod > p_enum->log_p_min) {
				p_enum->curr_entry = curr_entry + 1;
				return curr_prod->prod;
			}
		}

		if (p_enum->next_prime == s->aprog_data.num_aprogs)
			break;

		enum_run(s);
		p_enum->curr_entry = 0;
	}

	return P_SEARCH_DONE;
}

/*------------------------------------------------------------------------*/
uint64
sieve_fb_next(sieve_fb_t *s, poly_search_t *poly,
		root_callback callback, void *extra)
{
	uint32 i, j;
	uint64 p;
	uint32 num_roots;
	uint32 num_factors;
	uint32 factors[MAX_FACTORS];
	uint32 found;

	while (1) {

		if (s->curr_algo == ALGO_ENUM || 
		    s->curr_algo == ALGO_SIEVE) {

			if (s->curr_algo == ALGO_ENUM)
				p = get_next_enum_composite(s);
			else
				p = get_next_composite(s);

			if (p == P_SEARCH_DONE) {
				uint32 last_p = s->aprog_data.aprogs[
						s->aprog_data.num_aprogs-1].p;

				if (s->p_max >= P_PRIME_LIMIT ||
				    s->p_min >= P_PRIME_LIMIT ||
				    last_p >= s->p_max ||
				    s->num_roots_min > s->degree) {
					break;
				}

				s->curr_algo = ALGO_PRIME;
				init_prime_sieve(&s->p_prime, 
						MAX(s->p_min, last_p + 1),
						(uint32)s->p_max);
				continue;
			}

			num_factors = get_composite_factors(s, p, factors);

			for (i = found = 0; i < poly->num_poly; i++) {
				num_roots = get_composite_roots(s, 
							poly->batch + i, i, p, 
							num_factors, factors,
							s->num_roots_min,
							s->num_roots_max);

				if (num_roots == 0)
					continue;

				if (num_roots == INVALID_NUM_ROOTS)
					break;

				found++;
				lift_roots(s, poly->batch + i, p, num_roots);
				callback(p, num_roots, i, s->roots, extra);
			}

			if (found > 0)
				break;
		}
		else {
			p = get_next_prime(&s->p_prime);

			if (p >= s->p_max) {
				p = P_SEARCH_DONE;
				break;
			}

			for (i = found = 0; i < poly->num_poly; i++) {

				uint32 roots[MAX_POLYSELECT_DEGREE];

				num_roots = get_prime_roots(poly, i,
							(uint32)p, roots);

				if (num_roots == 0)
					continue;

				if (num_roots < s->num_roots_min ||
				    num_roots > s->num_roots_max)
					break;

				found++;
				for (j = 0; j < num_roots; j++) {
					mpz_set_ui(s->roots[j], 
						(mp_limb_t)roots[j]);
				}
				lift_roots(s, poly->batch + i, p, num_roots);
				callback(p, num_roots, i, s->roots, extra);
			}

			if (found > 0)
				break;
		}
	}

	return p;
}
