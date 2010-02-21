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

#include "stage1_core_deg5_96.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SHARED_BATCH_SIZE 32

typedef struct {
	uint64 p[SHARED_BATCH_SIZE];
	uint64 lattice_size[SHARED_BATCH_SIZE];
	uint32 roots[3 * POLY_BATCH_SIZE][SHARED_BATCH_SIZE];
} p_soa_shared_t;

__shared__ p_soa_shared_t pbatch_cache;

__constant__ uint96 two = {{2, 0, 0}};

/*------------------------------------------------------------------------*/
/* note that num_q must be a multiple of the block size
   (we want either all threads or no threads to execute
   the __syncthreads() call below) */

__global__ void
sieve_kernel(p_soa_t *pbatch, 
             uint32 num_p,
	     q_soa_t *qbatch,
	     uint32 num_q,
	     uint32 num_roots,
	     found_t *found_array)
{
	uint32 my_threadid;
	uint32 num_threads;
	uint32 i, j, k;

	my_threadid = blockIdx.x * blockDim.x + threadIdx.x;
	num_threads = gridDim.x * blockDim.x;
	found_array[my_threadid].p = 0;

	for (i = my_threadid; i < num_q; i += num_threads) {
		uint32 q2_w, p_done = 0;
		uint64 q;
		uint96 q2, q2_r;

		q = qbatch->p[i];
		q2 = wide_sqr48(q);
		q2_w = montmul32_w(q2.w[0]);
		q2_r = montmul96_r(q2, q2_w);

		while (p_done < num_p) {

			uint32 curr_num_p = MIN(SHARED_BATCH_SIZE,
						num_p - p_done);

			if (threadIdx.x < curr_num_p) {
				j = threadIdx.x;

				pbatch_cache.p[j] = pbatch->p[p_done + j];
				pbatch_cache.lattice_size[j] = 
					pbatch->lattice_size[p_done + j];

				for (k = 0; k < 3 * num_roots; k += 3) {
					pbatch_cache.roots[k][j] = 
						pbatch->roots[k][p_done + j];
					pbatch_cache.roots[k+1][j] = 
						pbatch->roots[k+1][p_done + j];
					pbatch_cache.roots[k+2][j] = 
						pbatch->roots[k+2][p_done + j];
				}
			}

			__syncthreads();

			for (j = 0; j < curr_num_p; j++) {
				uint32 prefetch0 = qbatch->roots[0][i];
				uint32 prefetch1 = qbatch->roots[1][i];
				uint32 prefetch2 = qbatch->roots[2][i];

				uint64 p = pbatch_cache.p[j];
				uint96 p2 = wide_sqr48(p);
				uint64 pinvmodq = modinv64(p, q);

				uint64 lattice_size = 
						pbatch_cache.lattice_size[j];
				uint96 pinv, tmp;
				uint96 test1;

				tmp = wide_sqr48(pinvmodq);
				tmp = montmul96(tmp, q2_r, q2, q2_w);
				pinv = montmul96(p2, tmp, q2, q2_w);
				pinv = modsub96(two, pinv, q2);
				pinv = montmul96(pinv, tmp, q2, q2_w);
				pinv = montmul96(pinv, q2_r, q2, q2_w);

				test1.w[0] = (uint32)lattice_size;
				test1.w[1] = (uint32)(lattice_size >> 32);
				test1.w[2] = 0;

				for (k = 0; k < 3 * num_roots; k += 3) {

					uint96 proot, qroot, res;

					qroot.w[0] = prefetch0;
					prefetch0 = qbatch->roots[k+3][i];
					qroot.w[1] = prefetch1;
					prefetch1 = qbatch->roots[k+4][i];
					qroot.w[2] = prefetch2;
					prefetch2 = qbatch->roots[k+5][i];

					proot.w[0] = pbatch_cache.roots[k][j];
					proot.w[1] = pbatch_cache.roots[k+1][j];
					proot.w[2] = pbatch_cache.roots[k+2][j];

					res = montmul96(pinv, 
							modsub96(qroot, proot,
							q2), q2, q2_w);

					if (cmp96(res, test1) < 0) {
						found_t *f = found_array + 
								my_threadid;
						f->p = p;
						f->q = q;
						f->which_poly = k / 3;
						f->offset = res;
						f->proot = proot;
					}
				}
			}

			p_done += curr_num_p;
		}
	}
}

#ifdef __cplusplus
}
#endif
