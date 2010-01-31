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
#include "stage1_core_deg46_64.h"

#define HOST_BATCH_SIZE 50000

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_roots;
	uint32 num_p;
	uint32 num_p_alloc;

	uint32 *p;
	uint64 *roots[MAX_ROOTS];
} q_soa_var_t;

#define MAX_P_SOA_ARRAYS 7

typedef struct {
	uint32 num_arrays;
	uint32 num_p;
	q_soa_var_t soa[MAX_P_SOA_ARRAYS];
} q_soa_array_t;

static void
q_soa_array_init(q_soa_array_t *s, uint32 degree)
{
	uint32 i, j;
	memset(s, 0, sizeof(q_soa_array_t));

	if (degree == 4) {
		s->num_arrays = 2;
		s->soa[0].num_roots = 8;
		s->soa[1].num_roots = 4;
	}
	else {  /* degree 6 */
		s->num_arrays = 7;
		s->soa[6].num_roots = 16;
		s->soa[5].num_roots = 24;
		s->soa[4].num_roots = 32;
		s->soa[3].num_roots = 36;
		s->soa[2].num_roots = 48;
		s->soa[1].num_roots = 64;
		s->soa[0].num_roots = 72;
	}

	for (i = 0; i < s->num_arrays; i++) {
		q_soa_var_t *soa = s->soa + i;

		soa->num_p_alloc = 1000;
		soa->p = (uint32 *)xmalloc(soa->num_p_alloc * 
					sizeof(uint32));
		for (j = 0; j < soa->num_roots; j++) {
			soa->roots[j] = (uint64 *)xmalloc(
						soa->num_p_alloc * 
						sizeof(uint64));
		}
	}
}

static void
q_soa_array_free(q_soa_array_t *s)
{
	uint32 i, j;

	for (i = 0; i < s->num_arrays; i++) {
		q_soa_var_t *soa = s->soa + i;

		free(soa->p);
		for (j = 0; j < soa->num_roots; j++)
			free(soa->roots[j]);
	}
}

static void
q_soa_array_reset(q_soa_array_t *s)
{
	uint32 i;

	s->num_p = 0;
	for (i = 0; i < s->num_arrays; i++)
		s->soa[i].num_p = 0;
}

static void
q_soa_var_grow(q_soa_var_t *soa)
{
	uint32 i;

	soa->num_p_alloc *= 2;
	soa->p = (uint32 *)xrealloc(soa->p, 
				soa->num_p_alloc * 
				sizeof(uint32));
	for (i = 0; i < soa->num_roots; i++) {
		soa->roots[i] = (uint64 *)xrealloc(soa->roots[i], 
					soa->num_p_alloc * 
					sizeof(uint64));
	}
}

static void 
store_p_soa(uint64 p, uint32 num_roots, uint32 which_poly,
		mpz_t *roots, void *extra)
{
	uint32 i, j;
	lattice_fb_t *L = (lattice_fb_t *)extra;
	q_soa_array_t *s = (q_soa_array_t *)(L->q_array);

	if (which_poly != 0) {
		printf("error: batched polynomials not allowed\n");
		exit(-1);
	}

	for (i = 0; i < s->num_arrays; i++) {
		uint32 num;
		q_soa_var_t *soa = s->soa + i;

		if (soa->num_roots != num_roots)
			continue;

		num = soa->num_p;
		if (soa->num_p_alloc == num)
			q_soa_var_grow(soa);

		soa->p[num] = (uint32)p;
		for (j = 0; j < num_roots; j++)
			soa->roots[j][num] = gmp2uint64(roots[j]);

		soa->num_p++;
		s->num_p++;
		break;
	}
}

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_p;
	uint32 p_size;
	uint32 p_size_alloc;
	p_packed_t *curr;
	p_packed_t *packed_array;
} p_packed_var_t;

static void 
p_packed_init(p_packed_var_t *s)
{
	memset(s, 0, sizeof(p_packed_var_t));

	s->p_size_alloc = 1000;
	s->packed_array = s->curr = (p_packed_t *)xmalloc(s->p_size_alloc *
						sizeof(p_packed_t));
}

static void 
p_packed_free(p_packed_var_t *s)
{
	free(s->packed_array);
}

static void 
p_packed_reset(p_packed_var_t *s)
{
	s->num_p = s->p_size = 0;
	s->curr = s->packed_array;
}

static p_packed_t * 
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + curr->num_roots);
}

static void 
store_p_packed(uint64 p, uint32 num_roots, uint32 which_poly,
		mpz_t *roots, void *extra)
{
	uint32 i;
	lattice_fb_t *L = (lattice_fb_t *)extra;
	p_packed_var_t *s = (p_packed_var_t *)(L->p_array);
	p_packed_t *curr;

	if (which_poly != 0) {
		printf("error: batched polynomials not allowed\n");
		exit(-1);
	}

	if ((p_packed_t *)((uint64 *)s->curr + s->p_size) + 1 >=
			s->packed_array + s->p_size_alloc ) {

		s->p_size_alloc *= 2;
		s->packed_array = (p_packed_t *)xrealloc(
						s->packed_array,
						s->p_size_alloc *
						sizeof(p_packed_t));
		s->curr = (p_packed_t *)((uint64 *)s->packed_array + s->p_size);
	}

	curr = s->curr;
	curr->p = (uint32)p;
	curr->lattice_size = 2 * L->poly->batch[0].sieve_size / 
				((double)p * p);
	curr->num_roots = num_roots;
	curr->pad = 0;
	for (i = 0; i < num_roots; i++)
		curr->roots[i] = gmp2uint64(roots[i]);

	s->num_p++;
	s->curr = p_packed_next(s->curr);
	s->p_size = ((uint8 *)s->curr - 
			(uint8 *)s->packed_array) / sizeof(uint64);
}

/*------------------------------------------------------------------------*/
static uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L,
			uint32 which_poly,
			uint32 threads_per_block,
			p_packed_var_t *p_array,
			q_soa_array_t *q_array,
			gpu_info_t *gpu_info, 
			CUfunction gpu_kernel,
			uint32 deadline)
{
	uint32 i, j;
	p_packed_t *packed_array = p_array->packed_array;
	q_soa_t *q_marshall = (q_soa_t *)L->q_marshall;
	uint32 num_blocks;
	uint32 num_p_offset;
	uint32 num_q_offset;
	uint32 num_qroots_offset;
	uint32 found_array_size = L->found_array_size;
	found_t *found_array = (found_t *)L->found_array;
	clock_t start_time = clock();
	void *gpu_ptr;

	i = 0;
	gpu_ptr = (void *)(size_t)L->gpu_q_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_q_offset = i;
	i += sizeof(uint32);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_qroots_offset = i;
	i += sizeof(uint32);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_p_offset = i;
	i += sizeof(uint32);

	gpu_ptr = (void *)(size_t)L->gpu_found_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_TRY(cuParamSetSize(gpu_kernel, i))


	for (i = 0; i < q_array->num_arrays; i++) {

		q_soa_var_t *soa = q_array->soa + i;
		uint32 num_qroots = soa->num_roots;
		uint32 num_q_done = 0;

		if (soa->num_p < 1000)
			continue;

		CUDA_TRY(cuParamSeti(gpu_kernel, 
				num_qroots_offset, num_qroots))

		while (num_q_done < soa->num_p) {

			uint32 num_p_done = 0;
			uint32 packed_words = 0;
			uint32 curr_num_p = 0;
			p_packed_t *packed_start = packed_array;

			uint32 curr_num_q = MIN(3 * found_array_size,
						soa->num_p - num_q_done);

			curr_num_q = MIN(curr_num_q, Q_SOA_BATCH_SIZE);

			memcpy(q_marshall->p, 
				soa->p + num_q_done,
				curr_num_q * sizeof(uint32));
			for (j = 0; j < num_qroots; j++) {
				memcpy(q_marshall->roots[j],
					soa->roots[j] + num_q_done,
					curr_num_q * sizeof(uint64));
			}

			CUDA_TRY(cuMemcpyHtoD(L->gpu_q_array, q_marshall,
					Q_SOA_BATCH_SIZE * (sizeof(uint32) +
						num_qroots * sizeof(uint64))))

			CUDA_TRY(cuParamSeti(gpu_kernel, num_q_offset, 
						curr_num_q))

			while (num_p_done < p_array->num_p) {
				p_packed_t *curr_packed = packed_start;

				do {
					uint32 next_words = packed_words +
							P_PACKED_HEADER_WORDS +
							curr_packed->num_roots;

					if (next_words >= P_ARRAY_WORDS)
						break;

					curr_num_p++;
					packed_words = next_words;
					curr_packed = p_packed_next(
								curr_packed);
				} while (++num_p_done < p_array->num_p);

#if 0
				printf("qroots %u qnum %u pnum %u pwords %u\n",
						num_qroots, curr_num_q,
						curr_num_p, packed_words);
#endif
				CUDA_TRY(cuMemcpyHtoD(L->gpu_p_array, 
							packed_start,
							packed_words *
							sizeof(uint64)))

				CUDA_TRY(cuParamSeti(gpu_kernel, num_p_offset, 
							curr_num_p))

				num_blocks = gpu_info->num_compute_units;
				if (curr_num_q < found_array_size) {
					num_blocks = (curr_num_q + 
						threads_per_block - 1) /
						threads_per_block;
				}

				CUDA_TRY(cuLaunchGrid(gpu_kernel, 
							num_blocks, 1))

				CUDA_TRY(cuMemcpyDtoH(found_array, 
							L->gpu_found_array, 
							threads_per_block * 
							num_blocks *
							sizeof(found_t)))

				for (j = 0; j < threads_per_block *
						num_blocks; j++) {
					found_t *f = found_array + j;

					if (f->p > 0) {
						uint128 proot = {{0}};
						uint128 offset = {{0}};

						proot.w[0] = (uint32)f->proot;
						proot.w[1] = 
						      (uint32)(f->proot >> 32);
						offset.w[0] = (uint32)f->offset;
						offset.w[1] = 
						      (uint32)(f->offset >> 32);

						handle_collision(L->poly, 
								which_poly,
								(uint64)f->p, 
								proot, offset, 
								(uint64)f->q);
					}
				}

				if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
					return 1;

				if ((double)(clock() - start_time) /
						CLOCKS_PER_SEC > deadline)
					return 0;

				packed_start = curr_packed;
				packed_words = 0;
				curr_num_p = 0;
			}

			num_q_done += curr_num_q;
		}
	}

	return 0;
}

/*------------------------------------------------------------------------*/
uint32
sieve_lattice_gpu_deg46_64(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel)
{
	uint32 i;
	uint32 min_small, min_large;
	uint32 quit = 0;
	p_packed_var_t * p_array;
	q_soa_array_t * q_array;
	uint32 degree = L->poly->degree;
	uint32 p_min_roots, p_max_roots;
	uint32 q_min_roots, q_max_roots;
	uint32 threads_per_block;

	L->q_marshall = (q_soa_t *)xmalloc(sizeof(q_soa_t));
	q_array = L->q_array = (q_soa_array_t *)xmalloc(
					sizeof(q_soa_array_t));
	p_array = L->p_array = (p_packed_var_t *)xmalloc(
					sizeof(p_packed_var_t));
	p_packed_init(p_array);
	q_soa_array_init(q_array, degree);

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

	printf("------- %u-%u %u-%u\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	if (degree == 4) {
		p_min_roots = 4;
		p_max_roots = 8;
		q_min_roots = 4;
		q_max_roots = 128;
	}
	else {
		p_min_roots = 12;
		p_max_roots = 36;
		q_min_roots = 16;
		q_max_roots = 128;
	}

	min_large = large_p_min;
	sieve_fb_reset(sieve_small, (uint64)large_p_min, 
			(uint64)large_p_max, q_min_roots, 
			q_max_roots);

	while (min_large < large_p_max) {

		q_soa_array_reset(q_array);

		for (i = 0; i < HOST_BATCH_SIZE && 
				min_large < large_p_max; i++) {
			min_large = sieve_fb_next(sieve_small, L->poly,
						store_p_soa, L);
		}
		if (q_array->num_p == 0)
			break;

		min_small = small_p_min;
		sieve_fb_reset(sieve_large, 
				(uint64)small_p_min, (uint64)small_p_max,
				p_min_roots, p_max_roots);

		while (min_small <= small_p_max) {

			time_t curr_time;
			double elapsed;

			p_packed_reset(p_array);

			for (i = 0; i < HOST_BATCH_SIZE && 
					min_small < small_p_max; i++) {
				min_small = sieve_fb_next(sieve_large, L->poly,
							store_p_packed, L);
			}
			if (p_array->num_p == 0)
				break;

			if (sieve_lattice_batch(obj, L, 0,
					threads_per_block,
					p_array, q_array,
					gpu_info, gpu_kernel,
					L->deadline)) {
				quit = 1;
				goto finished;
			}

			curr_time = time(NULL);
			elapsed = curr_time - L->start_time;
			if (elapsed > L->deadline) {
				quit = 1;
				goto finished;
			}
		}
	}

finished:
	CUDA_TRY(cuMemFree(L->gpu_q_array))
	CUDA_TRY(cuMemFree(L->gpu_found_array))
	p_packed_free(p_array);
	q_soa_array_free(q_array);
	free(p_array);
	free(q_array);
	free(L->found_array);
	free(L->q_marshall);
	return quit;
}


