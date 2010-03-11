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
#include "stage1_core_deg5_96.h"

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_p;
	uint64 last_p;

	uint64 *p;
	uint64 *lattice_size;
	uint32 *roots[3*POLY_BATCH_SIZE];
} p_soa_var_t;

static void
p_soa_var_init(p_soa_var_t *soa, uint32 batch_size)
{
	uint32 i;

	memset(soa, 0, sizeof(soa));
	soa->p = (uint64 *)xmalloc(batch_size * sizeof(uint64));
	soa->lattice_size = (uint64 *)xmalloc(batch_size * sizeof(uint64));
	for (i = 0; i < 3 * POLY_BATCH_SIZE; i += 3) {
		soa->roots[i] = (uint32 *)xmalloc(batch_size * 
						sizeof(uint32));
		soa->roots[i+1] = (uint32 *)xmalloc(batch_size * 
						sizeof(uint32));
		soa->roots[i+2] = (uint32 *)xmalloc(batch_size * 
						sizeof(uint32));
	}
}

static void
p_soa_var_free(p_soa_var_t *soa)
{
	uint32 i;

	free(soa->p);
	free(soa->lattice_size);
	for (i = 0; i < 3 * POLY_BATCH_SIZE; i += 3) {
		free(soa->roots[i]);
		free(soa->roots[i+1]);
		free(soa->roots[i+2]);
	}
}

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
	uint96 tmp = {{0}};

	if (num_roots != 1) {
		printf("error: num_roots > 1\n");
		exit(-1);
	}

	if (L->fill_p)
		soa = (p_soa_var_t *)L->p_array;
	else
		soa = (p_soa_var_t *)L->q_array;

	num = soa->num_p;
	mpz_export(&tmp, NULL, -1, sizeof(uint32), 0, 0, roots[0]);

	if (p != soa->last_p) {
		soa->p[num] = p;
		soa->num_p++;
		soa->last_p = p;
		soa->roots[3*which_poly][num] = tmp.w[0];
		soa->roots[3*which_poly+1][num] = tmp.w[1];
		soa->roots[3*which_poly+2][num] = tmp.w[2];
	}
	else {
		soa->roots[3*which_poly][num-1] = tmp.w[0];
		soa->roots[3*which_poly+1][num-1] = tmp.w[1];
		soa->roots[3*which_poly+2][num-1] = tmp.w[2];
	}
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L,
			uint32 threads_per_block,
			gpu_info_t *gpu_info, CUfunction gpu_kernel)
{
	uint32 i;
	p_soa_var_t * p_array = (p_soa_var_t *)L->p_array;
	p_soa_var_t * q_array = (p_soa_var_t *)L->q_array;
	uint32 num_blocks;
	uint32 num_p_offset;
	uint32 num_q_offset;
	void *gpu_ptr;
	uint32 num_poly = L->poly->num_poly;

	p_soa_t *p_marshall = (p_soa_t *)L->p_marshall;
	q_soa_t *q_marshall = (q_soa_t *)L->q_marshall;
	found_t *found_array = (found_t *)L->found_array;
	uint32 found_array_size = L->found_array_size;
	uint32 num_q_done = 0;

	i = 0;
	gpu_ptr = (void *)(size_t)L->gpu_p_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_p_offset = i;
	i += sizeof(uint32);

	gpu_ptr = (void *)(size_t)L->gpu_q_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_q_offset = i;
	i += sizeof(uint32);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	CUDA_TRY(cuParamSeti(gpu_kernel, i, num_poly))
	i += sizeof(uint32);

	gpu_ptr = (void *)(size_t)L->gpu_found_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_TRY(cuParamSetSize(gpu_kernel, i))

	while (num_q_done < q_array->num_p) {

		uint32 num_p_done = 0;
		time_t curr_time;
		double elapsed;
		uint32 curr_num_q = MIN(3 * found_array_size,
					q_array->num_p - num_q_done);

		curr_num_q = MIN(curr_num_q, Q_SOA_BATCH_SIZE);

		/* force to be a multiple of the block size */
		curr_num_q -= (curr_num_q % threads_per_block);
		if (curr_num_q == 0)
			break;

		memcpy(q_marshall->p, 
			q_array->p + num_q_done,
			curr_num_q * sizeof(uint64));

		for (i = 0; i < 3 * num_poly; i++) {
			memcpy(q_marshall->roots[i],
				q_array->roots[i] + num_q_done,
				curr_num_q * sizeof(uint32));
		}

		CUDA_TRY(cuMemcpyHtoD(L->gpu_q_array, q_marshall,
				Q_SOA_BATCH_SIZE * (sizeof(uint64) +
					3 * num_poly * sizeof(uint32))))
		CUDA_TRY(cuParamSeti(gpu_kernel, num_q_offset, curr_num_q))

		while (num_p_done < p_array->num_p) {

			uint32 curr_num_p = MIN(found_array_size / 3,
						p_array->num_p - num_p_done);

			curr_num_p = MIN(curr_num_p, P_SOA_BATCH_SIZE);
			memcpy(p_marshall->p, 
				p_array->p + num_p_done,
				curr_num_p * sizeof(uint64));
			memcpy(p_marshall->lattice_size, 
				p_array->lattice_size + num_p_done,
				curr_num_p * sizeof(uint64));

			for (i = 0; i < 3 * num_poly; i++) {
				memcpy(p_marshall->roots[i],
					p_array->roots[i] + num_p_done,
					curr_num_p * sizeof(uint32));
			}

			CUDA_TRY(cuMemcpyHtoD(L->gpu_p_array, p_marshall,
				P_SOA_BATCH_SIZE * (2 * sizeof(uint64) +
					3 * num_poly * sizeof(uint32))))

			CUDA_TRY(cuParamSeti(gpu_kernel, num_p_offset, 
						curr_num_p))
#if 0
			printf("qnum %u pnum %u\n", curr_num_q, curr_num_p);
#endif

			num_blocks = gpu_info->num_compute_units;
			if (curr_num_q < found_array_size)
				num_blocks = curr_num_q / threads_per_block;

			CUDA_TRY(cuLaunchGrid(gpu_kernel, 
						num_blocks, 1))

			CUDA_TRY(cuMemcpyDtoH(found_array, 
						L->gpu_found_array, 
						num_blocks * 
						threads_per_block *
							sizeof(found_t)))

			for (i = 0; i < threads_per_block * num_blocks; i++) {
				found_t *f = found_array + i;

				if (f->p > 0) {
					uint128 proot, offset;

					proot.w[0] = f->proot.w[0];
					proot.w[1] = f->proot.w[1];
					proot.w[2] = f->proot.w[2];
					proot.w[3] = 0;
					offset.w[0] = f->offset.w[0];
					offset.w[1] = f->offset.w[1];
					offset.w[2] = f->offset.w[2];
					offset.w[3] = 0;

					handle_collision(L->poly, 
						f->which_poly, f->p, proot, 
						offset, f->q);
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

		num_q_done += curr_num_q;
	}

	return 0;
}

/*------------------------------------------------------------------------*/
uint32
sieve_lattice_deg5_96(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max)
{
	uint32 i;
	uint64 min_small, min_large;
	uint32 quit = 0;
	p_soa_var_t * p_array;
	p_soa_var_t * q_array;
	uint32 num_poly = L->poly->num_poly;
	uint32 host_p_batch_size;
	uint32 host_q_batch_size;

	uint32 threads_per_block;
	gpu_info_t *gpu_info = L->gpu_info;
       	CUfunction gpu_kernel = L->gpu_kernel;

	L->p_marshall = (p_soa_t *)xmalloc(sizeof(p_soa_t));
	L->q_marshall = (q_soa_t *)xmalloc(sizeof(q_soa_t));
	p_array = L->p_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));
	q_array = L->q_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));

	CUDA_TRY(cuMemAlloc(&L->gpu_p_array, sizeof(p_soa_t)))
	CUDA_TRY(cuMemAlloc(&L->gpu_q_array, sizeof(q_soa_t)))

	CUDA_TRY(cuFuncGetAttribute((int *)&threads_per_block, 
			CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK,
			gpu_kernel))

	CUDA_TRY(cuFuncSetBlockShape(gpu_kernel, 
				threads_per_block, 1, 1))

	L->found_array_size = threads_per_block *
				gpu_info->num_compute_units;
	L->found_array = (found_t *)xmalloc(L->found_array_size *
					sizeof(found_t));
	CUDA_TRY(cuMemAlloc(&L->gpu_found_array, 
			L->found_array_size * sizeof(found_t)))

	host_p_batch_size = MAX(10000, L->found_array_size / 3);
	host_q_batch_size = MAX(50000, 12 * L->found_array_size);
	p_soa_var_init(p_array, host_p_batch_size);
	p_soa_var_init(q_array, host_q_batch_size);

	printf("------- %" PRIu64 "-%" PRIu64 " %" PRIu64 "-%" PRIu64 "\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	min_large = large_p_min;
	sieve_fb_reset(sieve_small, large_p_min, 
			large_p_max, 1, 1);

	while (min_large < large_p_max) {

		L->fill_p = 0;
		p_soa_var_reset(q_array);
		for (i = 0; i < host_q_batch_size && 
				min_large != P_SEARCH_DONE; i++) {
			min_large = sieve_fb_next(sieve_small, L->poly,
						store_p_soa, L);
		}
		if (q_array->num_p == 0)
			goto finished;

		min_small = small_p_min;
		sieve_fb_reset(sieve_large, 
				small_p_min, small_p_max,
				1, 1);

		while (min_small <= small_p_max) {

			L->fill_p = 1;
			p_soa_var_reset(p_array);
			for (i = 0; i < host_p_batch_size && 
					min_small != P_SEARCH_DONE; i++) {
				min_small = sieve_fb_next(sieve_large, L->poly,
							store_p_soa, L);
			}
			if (p_array->num_p == 0)
				goto finished;

			for (i = 0; i < p_array->num_p; i++) {
				uint64 p = p_array->p[i];

				p_array->lattice_size[i] = (uint64)
					(2 * L->poly->batch[num_poly - 
					1].sieve_size / ((double)p * p));
			}

			if (sieve_lattice_batch(obj, L, threads_per_block,
						gpu_info, gpu_kernel)) {
				quit = 1;
				goto finished;
			}
		}
	}

finished:
	CUDA_TRY(cuMemFree(L->gpu_p_array))
	CUDA_TRY(cuMemFree(L->gpu_q_array))
	CUDA_TRY(cuMemFree(L->gpu_found_array))
	p_soa_var_free(p_array);
	p_soa_var_free(q_array);
	free(p_array);
	free(q_array);
	free(L->p_marshall);
	free(L->q_marshall);
	free(L->found_array);
	return quit;
}


