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

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_p;
	uint32 last_p;

	uint32 p[HOST_BATCH_SIZE];
	uint32 lattice_size[HOST_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE][HOST_BATCH_SIZE];
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
	lattice_fb_t *L = (lattice_fb_t *)extra;
	p_soa_var_t *soa;
	uint32 num;

	if (num_roots != 1) {
		printf("error: num_roots > 1\n");
		exit(-1);
	}

	if (L->fill_p)
		soa = (p_soa_var_t *)L->p_array;
	else
		soa = (p_soa_var_t *)L->q_array;

	num = soa->num_p;
	if (p != soa->last_p) {
		soa->p[num] = (uint32)p;
		soa->num_p++;
		soa->last_p = (uint32)p;
		soa->roots[which_poly][num] = gmp2uint64(roots[0]);
	}
	else {
		soa->roots[which_poly][num - 1] = gmp2uint64(roots[0]);
	}
}

/*------------------------------------------------------------------------*/
static void
batch_invert(uint32 *plist, uint32 num_p, uint64 *invlist,
		uint64 q2, uint64 q2_r, uint32 q2_w)
{
	uint32 i;
	uint64 p2[INVERT_BATCH_SIZE];
	uint64 invprod;

	invlist[0] = invprod = wide_sqr32(plist[0]);
	for (i = 1; i < num_p; i++) {
		p2[i] = wide_sqr32(plist[i]);
		invlist[i] = invprod = montmul64(invprod, p2[i], q2, q2_w);
	}

	invprod = mp_modinv_2(invprod, q2);
	invprod = montmul64(invprod, q2_r, q2, q2_w);
	for (i = num_p - 1; i; i--) {
		invlist[i] = montmul64(invprod, invlist[i-1], q2, q2_w);
		invprod = montmul64(invprod, p2[i], q2, q2_w);
	}
	invlist[i] = invprod;
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L)
{
	uint32 i, j, k;
	p_soa_var_t * p_array = (p_soa_var_t *)L->p_array;
	p_soa_var_t * q_array = (p_soa_var_t *)L->q_array;
	uint32 num_poly = L->poly->num_poly;
	uint32 num_p = p_array->num_p;

	for (i = 0; i < q_array->num_p; i++) {

		uint32 q = q_array->p[i];
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
					q2, q2_r, q2_w);

			for (j = 0; j < curr_num_p; j++) {

				uint32 p = plist[j];
				uint64 pinv = pinvlist[j];
				uint32 lattice_size = 
					   p_array->lattice_size[num_p_done+j];

				for (k = 0; k < num_poly; k++) {

					uint64 qroot, proot, res;

					qroot = q_array->roots[k][i];
					proot = p_array->roots[k][num_p_done+j];
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

						handle_collision(L->poly, k,
								(uint64)p, 
								r, off, 
								(uint64)q);
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
sieve_lattice_deg5_64(msieve_obj *obj, lattice_fb_t *L, 
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
			(uint64)large_p_max, 1, 1);

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
				1, 1);

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


