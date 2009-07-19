#include "stage1.h"

#if 0
#define CHECK
#endif

typedef struct {
	p_batch_t p_array;
	p_batch_t q_array;

	poly_batch_t *poly;

	double start_time;
	uint32 deadline;
	uint32 num_tests;

	mpz_t p, p2, pr, q, q2, qr, inv, tmp, limit;
	mpz_t tmp1, tmp2, tmp3, tmp4;
} lattice_fb_t;

#define MIN_LATTICE_SIZE 8

/*------------------------------------------------------------------------*/
static void 
lattice_fb_free(lattice_fb_t *L)
{
	p_batch_free(&L->p_array);
	p_batch_free(&L->q_array);
	mpz_clear(L->p);
	mpz_clear(L->p2);
	mpz_clear(L->pr);
	mpz_clear(L->q);
	mpz_clear(L->q2);
	mpz_clear(L->qr);
	mpz_clear(L->inv);
	mpz_clear(L->tmp);
	mpz_clear(L->limit);
	mpz_clear(L->tmp1);
	mpz_clear(L->tmp2);
	mpz_clear(L->tmp3);
	mpz_clear(L->tmp4);
}

/*------------------------------------------------------------------------*/
static void 
lattice_fb_init(lattice_fb_t *L, poly_batch_t *poly, 
		uint32 deadline)
{
	p_batch_init(&L->p_array);
	p_batch_init(&L->q_array);
	L->poly = poly;

	L->start_time = get_cpu_time();
	L->deadline = deadline;
	L->num_tests = 0;

	mpz_init(L->p);
	mpz_init(L->p2);
	mpz_init(L->pr);
	mpz_init(L->q);
	mpz_init(L->q2);
	mpz_init(L->qr);
	mpz_init(L->inv);
	mpz_init(L->tmp);
	mpz_init(L->limit);
	mpz_init(L->tmp1);
	mpz_init(L->tmp2);
	mpz_init(L->tmp3);
	mpz_init(L->tmp4);
}

/*------------------------------------------------------------------------*/
static void
handle_collision(poly_batch_t *poly, uint32 which_poly,
		uint32 p, uint32 q, mpz_t m0_inc)
{
	curr_poly_t *curr = poly->batch + which_poly;

	mpz_set_ui(poly->p, (mp_limb_t)p);
	mpz_mul_ui(poly->p, poly->p, (mp_limb_t)q);
	mpz_add(poly->m0, curr->trans_m0, m0_inc);

#ifdef CHECK
	gmp_printf("poly %u p %u q %u coeff %Zd\n", which_poly, p, q, poly->p);

	mpz_pow_ui(poly->tmp1, poly->m0, (mp_limb_t)poly->degree);
	mpz_mul(poly->tmp2, poly->p, poly->p);
	mpz_sub(poly->tmp1, curr->trans_N, poly->tmp1);
	mpz_tdiv_r(poly->tmp3, poly->tmp1, poly->tmp2);
	if (mpz_cmp_ui(poly->tmp3, (mp_limb_t)0)) {
		printf("crap\n");
		exit(-1);
	}
#endif

	mpz_mul_ui(poly->tmp1, curr->high_coeff, (mp_limb_t)poly->degree);
	mpz_tdiv_qr(poly->m0, poly->tmp2, poly->m0, poly->tmp1);
	mpz_invert(poly->tmp3, poly->tmp1, poly->p);

	mpz_sub(poly->tmp4, poly->tmp3, poly->p);
	if (mpz_cmpabs(poly->tmp3, poly->tmp4) < 0)
		mpz_set(poly->tmp4, poly->tmp3);

	mpz_sub(poly->tmp5, poly->tmp2, poly->tmp1);
	if (mpz_cmpabs(poly->tmp2, poly->tmp5) > 0)
		mpz_add_ui(poly->m0, poly->m0, (mp_limb_t)1);
	else
		mpz_set(poly->tmp5, poly->tmp2);

	mpz_addmul(poly->m0, poly->tmp4, poly->tmp5);

	poly->callback(curr->high_coeff, poly->p, poly->m0, 
			curr->coeff_max, poly->callback_data);
}

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 lattice_size;
	uint64 *roots;
} cached_root_t;

static void
sieve_lattice_p32_q64(uint32 p, uint32 num_ppoly, 
		index_t *p_index, uint32 num_proots,
		uint64 *proots, uint32 q, uint32 num_qroots,
		cached_root_t cached_qroots[POLY_BATCH_SIZE],
		lattice_fb_t *L)
{
	uint32 i, j, k;
	poly_batch_t *poly = L->poly;
	uint32 p2 = p * p;
	uint64 q2 = (uint64)q * (uint64)q;
	uint32 inv = mp_modinv_1(mp_mod64(q2, p2), p2);

	for (i = 0; i < num_ppoly; i++) {

		uint32 which_poly = p_index[i].which_poly;
		cached_root_t *cache = cached_qroots + which_poly;
		uint64 *pr;
		uint64 *qr = cache->roots;
		uint32 lattice_size;

		if (qr == NULL)
			continue;

		pr = proots + p_index[i].start_offset;
		lattice_size = cache->lattice_size;

		for (j = 0; j < num_qroots; j++) {

			uint32 qrj = mp_mod64(qr[j], p2);

			for (k = 0; k < num_proots; k++) {
				uint32 prk = pr[k];
				uint32 tmp = mp_modmul_1(inv, 
						mp_modsub_1(prk, qrj, p2), p2);

				if (tmp < lattice_size ||
				    tmp >= p2 - lattice_size) {

					mpz_set_ui(L->tmp, (mp_limb_t)tmp);
					uint64_2gmp(q2, L->q2);
					uint64_2gmp(qr[j], L->qr);

					if (tmp > p2 - lattice_size) {
						mpz_sub_ui(L->tmp, L->tmp, 
							     (mp_limb_t)p2);
					}
					mpz_addmul(L->qr, L->tmp, L->q2);
					handle_collision(poly, which_poly,
						       	p, q, L->qr);
				}
			}
		}
	}
}

/*------------------------------------------------------------------------*/
static void
sieve_lattice_p64_q64(uint32 p, uint32 num_ppoly, 
		index_t *p_index, uint32 num_proots, 
		uint64 *proots, uint32 q, uint32 num_qroots,
		cached_root_t cached_qroots[POLY_BATCH_SIZE],
		lattice_fb_t *L)
{
	uint32 i, j, k;
	poly_batch_t *poly = L->poly;
	uint64 p2 = (uint64)p * (uint64)p;
	uint64 q2 = (uint64)q * (uint64)q;
	uint64 inv = mp_modinv_2(q2 % p2, p2);

#if !defined(GCC_ASM64A)
	uint64_2gmp(inv, L->inv);
	uint64_2gmp(p2, L->p2);
	uint64_2gmp(q2, L->q2);
#endif

	for (i = 0; i < num_ppoly; i++) {

		uint32 which_poly = p_index[i].which_poly;
		cached_root_t *cache = cached_qroots + which_poly;
		uint64 *pr;
		uint64 *qr = cache->roots;
		uint32 lattice_size;

		if (qr == NULL)
			continue;

		pr = proots + p_index[i].start_offset;
		lattice_size = cache->lattice_size;

		for (j = 0; j < num_qroots; j++) {

			uint64 qrj = qr[j] % p2;

			for (k = 0; k < num_proots; k++) {
				uint64 tmp = mp_modsub_2(pr[k], qrj, p2);

#if defined(GCC_ASM64A)		/*-----------------------------------*/
				ASM_G("mulq %1        \n\t"
				      "divq %2        \n\t"
				      "movq %%rdx, %0 \n\t"
				      : "+a"(tmp)
				      : "g"(inv), "g"(p2)
				      : "%rdx", "cc");

				if (tmp < lattice_size ||
				    tmp >= p2 - lattice_size) {

					uint64_2gmp(qr[j], L->qr);
					uint64_2gmp(tmp, L->tmp);

					if (tmp > p2 - lattice_size) {
						mpz_sub_ui(L->tmp, L->tmp, 
								(mp_limb_t)p2);
					}
					mpz_addmul_ui(L->qr, L->tmp, 
							(mp_limb_t)q2);
					handle_collision(poly, which_poly,
						       	p, q, L->qr);
				}
#else				/*-----------------------------------*/

	#if defined(_MSC_VER) && defined(_WIN64)
				uint64 mul_mod_64(uint64, uint64, uint64);

				tmp = mul_mod_64(tmp, inv, p2);
	#else
				uint64_2gmp(tmp, L->tmp);
				mpz_mul(L->tmp, L->tmp, L->inv);
				mpz_tdiv_r(L->tmp, L->tmp, L->p2);
				tmp = gmp2uint64(L->tmp);
	#endif

				if (tmp < lattice_size ||
				    tmp >= p2 - lattice_size) {

	#if defined(_MSC_VER) && defined(_WIN64)
					uint64_2gmp(tmp, L->tmp);
	#endif
					if (tmp > p2 - lattice_size) {
						mpz_sub(L->tmp, L->tmp, L->p2);
					}
					uint64_2gmp(qr[j], L->qr);
					mpz_addmul(L->qr, L->tmp, L->q2);
					handle_collision(poly, which_poly,
						       	p, q, L->qr);
				}
#endif				/*-----------------------------------*/
			}
		}
	}
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_q(lattice_fb_t *L)
{
	uint32 i;
	poly_batch_t *poly = L->poly;
	p_batch_t *p_array = &L->p_array;
	p_batch_t *q_array = &L->q_array;
	uint32 p_cutoff;

	p_entry_t *q_entry = q_array->entries;
	uint32 q = q_entry->p;
	uint32 num_qpoly = q_entry->num_poly;
	uint32 num_qroots = q_entry->num_roots;
	index_t *q_index = q_array->index;
	cached_root_t cached_qroots[POLY_BATCH_SIZE];
	uint32 max_lattice = 0;

	memset(cached_qroots, 0, sizeof(cached_qroots));

	for (i = 0; i < num_qpoly; i++) {
		uint32 which_poly = q_index[i].which_poly;
		curr_poly_t *curr = poly->batch + which_poly;
		uint32 lattice_size = (uint32)(curr->sieve_size / 
						((double)q * q));

		max_lattice = MAX(max_lattice, lattice_size);
		cached_qroots[which_poly].lattice_size = MAX(lattice_size, 2);
		cached_qroots[which_poly].roots = q_array->roots + 
						q_index[i].start_offset;
	}

	p_cutoff = (uint32)sqrt((double)((uint64)q * q >> 32)) + 1;

	for (i = 0; i < p_array->num_small_entries; i++) {

		p_entry_t *p_entry = p_array->entries + i;
		uint32 p = p_entry->p;
		uint32 num_ppoly = p_entry->num_poly;
		uint32 num_proots = p_entry->num_roots;
		index_t *p_index = p_array->index + p_entry->index_start_offset;
		uint64 *proots = p_array->roots;

		if (p > p_cutoff)
			sieve_lattice_p32_q64(p, num_ppoly, p_index, 
						num_proots, proots, 
						q, num_qroots, 
						cached_qroots, L);
		else
			sieve_lattice_p64_q64(p, num_ppoly, p_index, 
						num_proots, proots, 
						q, num_qroots, 
						cached_qroots, L);
	}

	for (; i < p_array->num_entries; i++) {

		p_entry_t *p_entry = p_array->entries + i;
		uint32 p = p_entry->p;
		uint32 num_ppoly = p_entry->num_poly;
		index_t *p_index = p_array->index + p_entry->index_start_offset;
		uint32 num_proots = p_entry->num_roots;
		uint64 *proots = p_array->roots;

		sieve_lattice_p64_q64(p, num_ppoly, p_index, 
					num_proots, proots,
					q, num_qroots, 
					cached_qroots, L);
	}

	L->num_tests += p_array->num_entries;
	if (L->num_tests >= 2000000) {

		double curr_time = get_cpu_time();
		double elapsed = curr_time - L->start_time;

		if (elapsed > L->deadline)
			return 1;
		L->num_tests = 0;
	}

	return (max_lattice < MIN_LATTICE_SIZE);
}

/*------------------------------------------------------------------------*/
#define P_BATCH_SIZE 5000

static uint32
sieve_lattice_core(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max)
{
	uint32 i;
	uint32 min_small, min_large;
	uint32 min_roots_small, min_roots_large;
	uint32 degree = L->poly->degree;

	printf("p %u %u %u %u lattice %.0lf\n", 
			small_p_min, small_p_max,
			large_p_min, large_p_max,
			L->poly->batch[L->poly->num_poly - 1].sieve_size / 
				((double)large_p_min * large_p_min));

	min_small = small_p_min;
	sieve_fb_reset(sieve_small, (uint64)min_small);

	if (small_p_min < 1000000)
		min_roots_small = 1;
	else if (small_p_min < 10000000)
		min_roots_small = degree;
	else
		min_roots_small = degree * degree;

	if (large_p_min < 1000000)
		min_roots_large = 1;
	else if (large_p_min < 3000000)
		min_roots_large = degree;
	else
		min_roots_large = degree * degree;

	if (degree == 4)
		min_roots_large *= 2;

	while (min_small < small_p_max) {

		p_batch_reset(&L->p_array);
		for (i = 0; i < P_BATCH_SIZE && 
					min_small < small_p_max; i++) {
			min_small = sieve_fb_next(sieve_small, L->poly,
						&L->p_array, 
						(uint64)small_p_max,
						min_roots_small);
		}
		if (L->p_array.num_entries == 0)
			return 0;

		printf("batch %u %u\n", L->p_array.num_entries, min_small);

		min_large = large_p_min;
		sieve_fb_reset(sieve_large, (uint64)min_large);

		while (min_large <= large_p_max) {

			p_batch_reset(&L->q_array);
			min_large = sieve_fb_next(sieve_large, L->poly,
						&L->q_array, 
						(uint64)large_p_max,
						min_roots_large);
			if (sieve_lattice_q(L) ||
			    obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				return 1;
		}
	}

	return 0;
}

/*------------------------------------------------------------------------*/
#define P_SCALE 1.3
#define MIN_SMALL_RANGE 20000

void
sieve_lattice(msieve_obj *obj, poly_batch_t *poly, 
		uint32 small_fb_max, uint32 large_fb_min, 
		uint32 large_fb_max, uint32 deadline)
{
	lattice_fb_t L;
	sieve_fb_t sieve_small, sieve_large;
	uint32 small_p_min, small_p_max;
	uint32 large_p_min, large_p_max;
	curr_poly_t *mid_poly = poly->batch + poly->num_poly / 2;
	curr_poly_t *last_poly = poly->batch + poly->num_poly - 1;

	if (last_poly->p_size_max >= (double)MAX_P * MAX_P) {
		printf("error: rational leading coefficient is too large\n");
		exit(-1);
	}
	large_p_min = sqrt(mid_poly->p_size_max);
	if (large_p_min >= MAX_P / P_SCALE)
		large_p_max = MAX_P - 1;
	else
		large_p_max = P_SCALE * large_p_min;

	small_p_min = large_p_min / P_SCALE;
	small_p_max = large_p_min - 1;

	gmp_printf("coeff %Zd-%Zd %u %u %u %u lattice %.0lf\n", 
			poly->batch[0].high_coeff, 
			poly->batch[poly->num_poly - 1].high_coeff,
			small_p_min, small_p_max,
			large_p_min, large_p_max,
			last_poly->sieve_size / 
				((double)large_p_min * large_p_min));

	if (small_p_max < large_fb_min * large_fb_min) {
		if (small_p_min < large_fb_max) {
			small_p_max = MIN(small_p_max, large_fb_max);
		}
		else {
			small_p_max = large_fb_max;
			while (small_p_min > large_fb_max)
				small_p_min = small_p_min / P_SCALE;
		}

		while (small_p_min > large_fb_min &&
		       small_p_max - small_p_min < MIN_SMALL_RANGE) {
			small_p_min = small_p_min / P_SCALE;
		}
		small_p_min = MAX(small_p_min, large_fb_min);
		large_p_min = mid_poly->p_size_max / small_p_max;
		large_p_max = mid_poly->p_size_max / small_p_min;
	}
	else if (small_p_min > large_fb_max &&
		 small_p_min < large_fb_min * large_fb_min) {
		small_p_min = large_fb_min * large_fb_min;
	}

	if (last_poly->sieve_size / ((double)large_p_min * 
			large_p_min) < MIN_LATTICE_SIZE) {
		return;
	}

	sieve_fb_init(&sieve_small, poly, large_fb_min, large_fb_max);
	sieve_fb_init(&sieve_large, poly, 5, small_fb_max);
	lattice_fb_init(&L, poly, deadline);

	while (1) {
		if (sieve_lattice_core(obj, &L, &sieve_small, &sieve_large,
					small_p_min, small_p_max,
					large_p_min, large_p_max)) {
			break;
		}

		small_p_max = small_p_min;
		if (small_p_max <= large_fb_min)
			break;
		while (small_p_min > large_fb_min &&
		       small_p_max - small_p_min < MIN_SMALL_RANGE) {
			small_p_min = small_p_min / P_SCALE;
		}
		small_p_max--;
		small_p_min = MAX(small_p_min, large_fb_min);

		if (small_p_max < large_fb_min * large_fb_min) {
			if (small_p_min < large_fb_max) {
				small_p_max = MIN(small_p_max, large_fb_max);
			}
			else {
				small_p_max = large_fb_max;
				while (small_p_min > large_fb_max)
					small_p_min = small_p_min / P_SCALE;
			}

			while (small_p_min > large_fb_min &&
			       small_p_max - small_p_min < MIN_SMALL_RANGE) {
				small_p_min = small_p_min / P_SCALE;
			}
			small_p_min = MAX(small_p_min, large_fb_min);
		}
		else if (small_p_min > large_fb_max &&
			 small_p_min < large_fb_min * large_fb_min) {
			small_p_min = large_fb_min * large_fb_min;
		}

		if (mid_poly->p_size_max / small_p_max >= MAX_P)
			break;
		large_p_min = mid_poly->p_size_max / small_p_max;

		if (mid_poly->p_size_max / small_p_min >= MAX_P)
			large_p_max = MAX_P - 1;
		else
			large_p_max = mid_poly->p_size_max / small_p_min;
	}

	lattice_fb_free(&L);
	sieve_fb_free(&sieve_small);
	sieve_fb_free(&sieve_large);
}
