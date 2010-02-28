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

	uint64 p[HOST_BATCH_SIZE];
	uint32 num_roots[HOST_BATCH_SIZE];
	uint64 lattice_size[HOST_BATCH_SIZE];
	uint128 roots[MAX_ROOTS][HOST_BATCH_SIZE];
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
	soa->p[num] = p;
	soa->num_roots[num] = num_roots;
	for (i = 0; i < num_roots; i++) {
		uint128 tmp = {{0}};

		mpz_export(&tmp, NULL, -1, sizeof(uint32), 0, 0, roots[i]);
		soa->roots[i][num] = tmp;
	}

	soa->num_p++;
}

static void
batch_invert(uint64 *plist, uint32 num_p, uint128 *invlist,
		uint128 q2, uint128 q2_r, uint32 q2_w)
{
	uint32 i;
	uint128 p2[INVERT_BATCH_SIZE];
	uint128 invprod;

	invlist[0] = invprod = wide_sqr64(plist[0]);
	for (i = 1; i < num_p; i++) {
		p2[i] = wide_sqr64(plist[i]);
		invlist[i] = invprod = montmul128(invprod, p2[i], q2, q2_w);
	}

	{
		mp_t mp_prod, mp_q2, res;

		mp_q2.nwords = 4;
		mp_q2.val[0] = q2.w[0];
		mp_q2.val[1] = q2.w[1];
		mp_q2.val[2] = q2.w[2];
		mp_q2.val[3] = q2.w[3];
		if (mp_q2.val[3] == 0)
			mp_q2.nwords--;

		mp_prod.nwords = 4;
		mp_prod.val[0] = invprod.w[0];
		mp_prod.val[1] = invprod.w[1];
		mp_prod.val[2] = invprod.w[2];
		mp_prod.val[3] = invprod.w[3];
		if (mp_prod.val[3] == 0)
			mp_prod.nwords--;

		mp_modinv(&mp_prod, &mp_q2, &res);
		invprod.w[0] = res.val[0];
		invprod.w[1] = res.val[1];
		invprod.w[2] = res.val[2];
		invprod.w[3] = res.val[3];
		invprod = montmul128(invprod, q2_r, q2, q2_w);
	}

	for (i = num_p - 1; i; i--) {
		invlist[i] = montmul128(invprod, invlist[i-1], q2, q2_w);
		invprod = montmul128(invprod, p2[i], q2, q2_w);
	}
	invlist[i] = invprod;
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L)
{
	uint32 i, j, k, m;
	p_soa_var_t * p_array = (p_soa_var_t *)L->p_array;
	p_soa_var_t * q_array = (p_soa_var_t *)L->q_array;
	uint32 num_p = p_array->num_p;

	for (i = 0; i < q_array->num_p; i++) {

		uint64 q = q_array->p[i];
		uint32 num_qroots = q_array->num_roots[i];
		uint128 q2 = wide_sqr64(q);
		uint32 q2_w = montmul32_w(q2.w[0]);
		uint128 q2_r = montmul128_r(q2);

		uint32 num_p_done = 0;
		time_t curr_time;
		double elapsed;

		while (num_p_done < num_p) {

			uint64 *plist = p_array->p + num_p_done;

			uint128 pinvlist[INVERT_BATCH_SIZE];
			uint32 curr_num_p = MIN(INVERT_BATCH_SIZE,
						num_p - num_p_done);

			batch_invert(plist, curr_num_p, pinvlist,
					q2, q2_r, q2_w);

			for (j = 0; j < curr_num_p; j++) {

				uint64 p = plist[j];
				uint128 pinv = pinvlist[j];
				uint32 num_proots = 
					  p_array->num_roots[num_p_done+j];
				uint64 lattice_size = 
					   p_array->lattice_size[num_p_done+j];
				uint128 test1;

				test1.w[0] = (uint32)lattice_size;
				test1.w[1] = (uint32)(lattice_size >> 32);
				test1.w[2] = 0;
				test1.w[3] = 0;

				for (k = 0; k < num_qroots; k++) {

				    uint128 qroot = q_array->roots[k][i];

				    for (m = 0; m < num_proots; m++) {

					uint128 proot, res;

					proot = p_array->roots[m][num_p_done+j];
					res = montmul128(pinv, 
							modsub128(qroot, proot,
							q2), q2, q2_w);

					if (cmp128(res, test1) < 0) {
						handle_collision(L->poly, 0,
								p, proot, 
								res, q);
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
sieve_lattice_deg6_128(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max)
{
	uint32 i;
	uint64 min_small, min_large;
	p_soa_var_t * p_array;
	p_soa_var_t * q_array;
	uint32 quit = 0;
	uint32 p_min_roots = 12;
	uint32 p_max_roots = 36;
	uint32 q_min_roots = 32;
	uint32 q_max_roots = 96;

	p_array = L->p_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));
	q_array = L->q_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));

	printf("------- %" PRIu64 "-%" PRIu64 " %" PRIu64 "-%" PRIu64 "\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	min_large = large_p_min;
	sieve_fb_reset(sieve_small, large_p_min, 
			large_p_max, q_min_roots, q_max_roots);

	while (min_large < large_p_max) {

		L->fill_p = 0;
		p_soa_var_reset(q_array);
		for (i = 0; i < HOST_BATCH_SIZE && 
				min_large != P_SEARCH_DONE; i++) {
			min_large = sieve_fb_next(sieve_small, L->poly,
						store_p_soa, L);
		}
		if (q_array->num_p == 0)
			goto finished;

		min_small = small_p_min;
		sieve_fb_reset(sieve_large, small_p_min, small_p_max,
				p_min_roots, p_max_roots);

		while (min_small <= small_p_max) {

			L->fill_p = 1;
			p_soa_var_reset(p_array);
			for (i = 0; i < HOST_BATCH_SIZE && 
					min_small != P_SEARCH_DONE; i++) {
				min_small = sieve_fb_next(sieve_large, L->poly,
							store_p_soa, L);
			}
			if (p_array->num_p == 0)
				goto finished;

			for (i = 0; i < p_array->num_p; i++) {
				uint64 p = p_array->p[i];

				p_array->lattice_size[i] = (uint64)
					(2 * L->poly->batch[0].sieve_size /  
					((double)p * p));
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
