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
#include "cpu_intrinsics.h"

#define HOST_BATCH_SIZE 10000

#define INVERT_BATCH_SIZE 200

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_p;
	uint32 last_p;

	uint32 p[HOST_BATCH_SIZE];
	uint32 num_roots[HOST_BATCH_SIZE];
	uint32 lattice_size[HOST_BATCH_SIZE];
	uint64 roots[MAX_ROOTS][HOST_BATCH_SIZE];
} p_soa_var_t;

static void
p_soa_var_reset(p_soa_var_t *soa)
{
	soa->num_p = 0;
	soa->last_p = 0;
}

static void 
store_p_soa(uint64 p, uint32 num_roots, uint32 which_poly, 
		mpz_t *roots, void *extra)
{
	uint32 i;
	lattice_fb_t *L = (lattice_fb_t *)extra;
	p_soa_var_t *soa;
	uint32 num;

	if (which_poly != 0) {
		printf("error: batched polynomials not allowed\n");
		exit(-1);
	}

	if (L->fill_p)
		soa = (p_soa_var_t *)L->p_array;
	else
		soa = (p_soa_var_t *)L->q_array;

	num = soa->num_p;
	soa->p[num] = (uint32)p;
	soa->num_roots[num] = num_roots;
	for (i = 0; i < num_roots; i++)
		soa->roots[i][num] = gmp2uint64(roots[i]);

	soa->num_p++;
}

/*------------------------------------------------------------------------*/
static void
batch_invert(uint32 *plist, uint32 num_p, uint64 *invlist,
		uint32 q, uint64 q2, uint64 q2_r, uint32 q2_w)
{
	uint32 i;

	for (i = 0; i < num_p; i++) {
		uint32 p = plist[i];
		uint64 p2 = wide_sqr32(p);
		uint64 inv = mp_modinv_2(p2, q2);
		invlist[i] = montmul64(inv, q2_r, q2, q2_w);
	}
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L)
{
	uint32 i, j, k, m;
	p_soa_var_t * p_array = (p_soa_var_t *)L->p_array;
	p_soa_var_t * q_array = (p_soa_var_t *)L->q_array;
	uint32 num_poly = L->poly->num_poly;
	uint32 num_p = p_array->num_p;

	for (i = 0; i < q_array->num_p; i++) {

		uint32 q = q_array->p[i];
		uint32 num_qroots = q_array->num_roots[i];
		uint64 q2 = wide_sqr32(q);
		uint32 q2_w = montmul32_w((uint32)q2);
		uint64 q2_r = montmul64_r(q2);

		uint32 num_p_done = 0;
		time_t curr_time;
		double elapsed;

		while (num_p_done < num_p) {

			uint32 *plist = p_array->p + num_p_done;

			uint64 pinvlist[INVERT_BATCH_SIZE];
			uint32 curr_num_p = MIN(INVERT_BATCH_SIZE,
						num_p - num_p_done);

			batch_invert(plist, curr_num_p, pinvlist,
					q, q2, q2_r, q2_w);

			for (j = 0; j < curr_num_p; j++) {

				uint32 p = plist[j];
				uint64 p2 = wide_sqr32(p);
				uint64 pinv = pinvlist[j];
				uint32 num_proots = 
					  p_array->num_roots[num_p_done+j];
				uint32 lattice_size = 
					   p_array->lattice_size[num_p_done+j];

				for (k = 0; k < num_qroots; k++) {

				    uint64 qroot = q_array->roots[k][i];

				    for (m = 0; m < num_proots; m++) {

					uint64 proot, res;

					proot = p_array->roots[m][num_p_done+j];
					res = montmul64(pinv, 
							modsub64(qroot, proot,
							q2), q2, q2_w);

					if (res < lattice_size) {
						uint128 r, off;

						r.w[0] = (uint32)proot;
						r.w[1] = (uint32)(proot >> 32);
						r.w[2] = 0;
						r.w[3] = 0;
						off.w[0] = (uint32)res;
						off.w[1] = (uint32)(res >> 32);
						off.w[2] = 0;
						off.w[3] = 0;

						handle_collision(L->poly, 0,
								(uint64)p, 
								r, off, 
								(uint64)q);
					}
				    }
				}
			}

			num_p_done += curr_num_p;
		}

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			return 1;

		curr_time = time(NULL);
		elapsed = curr_time - L->start_time;
		if (elapsed > L->deadline)
			return 1;
	}

	return 0;
}

/*------------------------------------------------------------------------*/
uint32
sieve_lattice_deg46_64(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max)
{
	uint32 i;
	uint32 min_small, min_large;
	uint32 quit = 0;
	p_soa_var_t * p_array;
	p_soa_var_t * q_array;
	uint32 num_poly = L->poly->num_poly;

	p_array = L->p_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));
	q_array = L->q_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));

	printf("------- %u-%u %u-%u\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	min_large = large_p_min;
	sieve_fb_reset(sieve_small, (uint64)large_p_min, 
			(uint64)large_p_max, 4, MAX_ROOTS);

	while (min_large < large_p_max) {

		L->fill_p = 0;
		p_soa_var_reset(q_array);
		for (i = 0; i < HOST_BATCH_SIZE && 
				min_large < large_p_max; i++) {
			min_large = sieve_fb_next(sieve_small, L->poly,
						store_p_soa, L);
		}
		if (q_array->num_p == 0)
			goto finished;

		min_small = small_p_min;
		sieve_fb_reset(sieve_large, 
				(uint64)small_p_min, (uint64)small_p_max,
				4, MAX_ROOTS);

		while (min_small <= small_p_max) {

			L->fill_p = 1;
			p_soa_var_reset(p_array);
			for (i = 0; i < HOST_BATCH_SIZE && 
					min_small < small_p_max; i++) {
				min_small = sieve_fb_next(sieve_large, L->poly,
							store_p_soa, L);
			}
			if (p_array->num_p == 0)
				goto finished;

			for (i = 0; i < p_array->num_p; i++) {
				uint32 p = p_array->p[i];

				p_array->lattice_size[i] = (uint32)
					(2 * L->poly->batch[num_poly - 
					1].sieve_size / ((double)p * p));
			}

			if (sieve_lattice_batch(obj, L)) {
				quit = 1;
				goto finished;
			}
		}
	}

finished:
	free(p_array);
	free(q_array);
	return quit;
}


