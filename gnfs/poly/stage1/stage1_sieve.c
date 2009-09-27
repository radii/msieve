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

#if 1
#define CHECK
#endif

/* structures for storing arithmetic progressions. Rational
   leading coeffs of NFS polynomials are assumed to be the 
   product of two groups of factors p, each of size at most 32
   bits (32 bits is enough for 512-bit factorizations), and 
   candidates must satisfy a condition modulo p^2 */

#define MAX_P_BITS 32
#define MAX_P ((uint64)1 << MAX_P_BITS)

#define P_BATCH_SIZE 1024

typedef struct {
	uint32 num_p;
	uint32 pad[15];

	uint32 num_roots[P_BATCH_SIZE];
	uint32 p[P_BATCH_SIZE];
	uint32 roots[2*MAX_ROOTS][P_BATCH_SIZE];
} p_batch_t;

typedef struct {
	p_batch_t *p_array;
	p_batch_t *q_array;

	poly_search_t *poly;

	double start_time;
	uint32 deadline;
	uint32 num_tests;
	uint32 fill_p_array;
} lattice_fb_t;

/*------------------------------------------------------------------------*/
static void 
lattice_fb_free(lattice_fb_t *L)
{
	aligned_free(L->p_array);
	aligned_free(L->q_array);
}

/*------------------------------------------------------------------------*/
static void 
lattice_fb_init(lattice_fb_t *L, poly_search_t *poly, 
		uint32 deadline)
{
	L->poly = poly;
	L->start_time = get_cpu_time();
	L->deadline = deadline;
	L->num_tests = 0;

	L->p_array = (p_batch_t *)aligned_malloc(sizeof(p_batch_t), 64);
	L->q_array = (p_batch_t *)aligned_malloc(sizeof(p_batch_t), 64);
}

/*------------------------------------------------------------------------*/
static void
handle_collision(poly_search_t *poly,
		uint64 p, uint64 proot, uint64 res, uint64 q)
{
	uint64_2gmp(p, poly->tmp1);
	uint64_2gmp(q, poly->tmp2);
	uint64_2gmp(proot, poly->tmp3);
	uint64_2gmp(res, poly->tmp4);

	mpz_mul(poly->p, poly->tmp1, poly->tmp2);
	mpz_mul(poly->tmp1, poly->tmp1, poly->tmp1);
	if (mpz_cmp(poly->tmp4, poly->tmp2) > 0)
		mpz_submul(poly->tmp4, poly->tmp2, poly->tmp2);

	mpz_addmul(poly->tmp3, poly->tmp1, poly->tmp4);
	mpz_add(poly->m0, poly->trans_m0, poly->tmp3);

#ifdef CHECK
	gmp_printf("p %.0lf q %.0lf coeff %Zd\n", 
				(double)p, (double)q, poly->p);

	mpz_pow_ui(poly->tmp1, poly->m0, (mp_limb_t)poly->degree);
	mpz_mul(poly->tmp2, poly->p, poly->p);
	mpz_sub(poly->tmp1, poly->trans_N, poly->tmp1);
	mpz_tdiv_r(poly->tmp3, poly->tmp1, poly->tmp2);
	if (mpz_cmp_ui(poly->tmp3, (mp_limb_t)0)) {
		printf("crap\n");
		exit(-1);
	}
#endif

	mpz_mul_ui(poly->tmp1, poly->high_coeff, (mp_limb_t)poly->degree);
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

	poly->callback(poly->high_coeff, poly->p, poly->m0, 
			poly->coeff_max, poly->callback_data);
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(lattice_fb_t *L)
{
	uint32 i, j, k, m;
	double sieve_size = L->poly->sieve_size;
	p_batch_t *pbatch = L->p_array;
	p_batch_t *qbatch = L->q_array;
	uint32 num_p = pbatch->num_p;
	uint32 num_q = qbatch->num_p;

	uint32 lattice_size[P_BATCH_SIZE];
	uint64 inv[P_BATCH_SIZE];

	for (i = 0; i < num_p; i++) {
		uint32 p = pbatch->p[i];
		lattice_size[i] = sieve_size / ((double)p * p);
	}

	for (i = 0; i < num_q; i++) {
		uint32 q = qbatch->p[i];
		uint64 q2 = (uint64)q * q;
		uint32 num_qroots = qbatch->num_roots[i];
		uint64 qroots[MAX_ROOTS];

		for (j = 0; j < num_qroots; j++) {
			uint32 qr0 = qbatch->roots[2*j][i];
			uint32 qr1 = qbatch->roots[2*j+1][i];
			qroots[j] = (uint64)qr1 << 32 | qr0;
		}

		for (j = 0; j < num_p; j++) {
			uint32 p = pbatch->p[j];
			uint64 p2 = (uint64)p * p;
			inv[j] = mp_modinv_2(p2, q2);
		}

		for (j = 0; j < num_p; j++) {
			uint32 num_proots = pbatch->num_roots[j];
			uint32 plattice = lattice_size[j];
			uint64 pinv = inv[j];

			for (k = 0; k < num_proots; k++) {

				uint32 pr0 = pbatch->roots[2*k][j];
				uint32 pr1 = pbatch->roots[2*k+1][j];
				uint64 proot = (uint64)pr1 << 32 | pr0;
						
				for (m = 0; m < num_qroots; m++) {
	
					uint64 res = mp_modmul_2(pinv,
								mp_modsub_2(
								   qroots[m],
								   proot, q2),
								q2);

					if (res < plattice ||
					    res >= q2 - plattice) {
						handle_collision(L->poly,
							(uint64)pbatch->p[j],
							proot, res, (uint64)q);
					}
				}
			}
		}
	}

	return 0;
}

/*------------------------------------------------------------------------*/
static void 
store_progression(uint64 p, uint32 num_roots,
		mpz_t *roots, void *extra)
{
	lattice_fb_t *L = (lattice_fb_t *)extra;
	p_batch_t *batch;
	uint32 num;
	uint32 i;

	batch = L->q_array;
	if (L->fill_p_array)
		batch = L->p_array;

	num = batch->num_p;
	batch->p[num] = (uint32)p;
	batch->num_roots[num] = num_roots;

	for (i = 0; i < num_roots; i++) {
		uint64 root = gmp2uint64(roots[i]);

		batch->roots[2*i][num] = (uint32)root;
		batch->roots[2*i+1][num] = (uint32)(root >> 32);
	}
	batch->num_p++;
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_core(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max)
{
	uint32 i;
	uint32 min_small, min_large;

	printf("p %u %u %u %u\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	min_small = small_p_min;
	sieve_fb_reset(sieve_small, 
			(uint64)small_p_min, (uint64)small_p_max,
			1, MAX_ROOTS);

	while (min_small < small_p_max) {

		L->fill_p_array = 1;
		L->p_array->num_p = 0;
		for (i = 0; i < P_BATCH_SIZE && 
				min_small < small_p_max; i++) {
			min_small = sieve_fb_next(sieve_small, L->poly,
						store_progression, L);
		}
		if (L->p_array->num_p == 0)
			return 0;

		printf("batch %u %u\n", L->p_array->num_p, min_small);

		min_large = large_p_min;
		sieve_fb_reset(sieve_large, 
				(uint64)large_p_min, (uint64)large_p_max,
				1, MAX_ROOTS);

		while (min_large <= large_p_max) {

			L->fill_p_array = 0;
			L->q_array->num_p = 0;
			for (i = 0; i < P_BATCH_SIZE && 
					min_large < large_p_max; i++) {
				min_large = sieve_fb_next(sieve_large, L->poly,
							store_progression, L);
			}
			if (L->q_array->num_p == 0)
				break;

			if (sieve_lattice_batch(L) ||
			    obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				return 1;
		}
	}

	return 0;
}

/*------------------------------------------------------------------------*/
#define P_SCALE 1.1

void
sieve_lattice(msieve_obj *obj, poly_search_t *poly, 
		uint32 small_fb_max, uint32 large_fb_min, 
		uint32 large_fb_max, uint32 deadline)
{
	lattice_fb_t L;
	sieve_fb_t sieve_small, sieve_large;
	uint32 small_p_min, small_p_max;
	uint32 large_p_min, large_p_max;

	if (poly->p_size_max >= (double)MAX_P * MAX_P) {
		printf("error: rational leading coefficient is too large\n");
		exit(-1);
	}
	large_p_min = sqrt(poly->p_size_max);
	if (large_p_min >= MAX_P / P_SCALE)
		large_p_max = MAX_P - 1;
	else
		large_p_max = P_SCALE * large_p_min;

	small_p_min = large_p_min / P_SCALE;
	small_p_max = large_p_min - 1;

	gmp_printf("coeff %Zd"
		   " %" PRIu64 " %" PRIu64 
		   " %" PRIu64 " %" PRIu64 "\n",
			poly->high_coeff, 
			(uint64)small_p_min, (uint64)small_p_max,
			(uint64)large_p_min, (uint64)large_p_max);

	sieve_fb_init(&sieve_small, poly, large_fb_min, large_fb_max);
	sieve_fb_init(&sieve_large, poly, 5, small_fb_max);
	lattice_fb_init(&L, poly, deadline);

	while (1) {
		if (sieve_lattice_core(obj, &L, &sieve_small, &sieve_large,
					small_p_min, small_p_max,
					large_p_min, large_p_max)) {
			break;
		}

		small_p_max = small_p_min - 1;
		small_p_min = small_p_min / P_SCALE;

		if (poly->p_size_max / small_p_max >= MAX_P)
			break;
		large_p_min = poly->p_size_max / small_p_max;

		if (poly->p_size_max / small_p_min >= MAX_P)
			large_p_max = MAX_P - 1;
		else
			large_p_max = poly->p_size_max / small_p_min;
	}

	lattice_fb_free(&L);
	sieve_fb_free(&sieve_small);
	sieve_fb_free(&sieve_large);
}
