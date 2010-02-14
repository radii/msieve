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

#include "stage1_core_deg5_64.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SHARED_BATCH_SIZE 48

typedef struct {
	uint32 p[SHARED_BATCH_SIZE];
	uint32 lattice_size[SHARED_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE][SHARED_BATCH_SIZE];
} p_soa_shared_t;

__shared__ p_soa_shared_t pbatch_cache;

/*------------------------------------------------------------------------*/
/* note that num_q must be a multiple of the block size
   (we want either all threads or no threads to execute
   the __syncthreads() call below) */

__global__ void
sieve_kernel_64(p_soa_t *pbatch, 
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
		uint32 q, q2_w, p_done = 0;
		uint64 q2, q2_r;

		q = qbatch->p[i];
		q2 = (uint64)q * q;
		q2_w = montmul32_w((uint32)q2);
		q2_r = montmul64_r(q2, q2_w);

		while (p_done < num_p) {

			uint32 curr_num_p = MIN(SHARED_BATCH_SIZE,
						num_p - p_done);

			if (threadIdx.x < curr_num_p) {
				j = threadIdx.x;

				pbatch_cache.p[j] = pbatch->p[p_done + j];
				pbatch_cache.lattice_size[j] = 
					pbatch->lattice_size[p_done + j];

				for (k = 0; k < num_roots; k++) {
					pbatch_cache.roots[k][j] = 
						pbatch->roots[k][p_done + j];
				}
			}

			__syncthreads();

			for (j = 0; j < curr_num_p; j++) {
				uint64 prefetch = qbatch->roots[0][i];
				uint32 p = pbatch_cache.p[j];
				uint64 p2 = (uint64)p * p;
				uint32 pinvmodq = modinv32(p, q);

				uint32 lattice_size = 
						pbatch_cache.lattice_size[j];
				uint64 pinv, tmp;

				tmp = (uint64)pinvmodq * pinvmodq;
				tmp = montmul64(tmp, q2_r, q2, q2_w);
				pinv = montmul64(p2, tmp, q2, q2_w);
				pinv = modsub64((uint64)2, pinv, q2);
				pinv = montmul64(pinv, tmp, q2, q2_w);
				pinv = montmul64(pinv, q2_r, q2, q2_w);

				for (k = 0; k < num_roots; k++) {

					uint64 qroot;
					uint64 proot;
					uint64 res;

					qroot = prefetch;
					prefetch = qbatch->roots[k+1][i];
					proot = pbatch_cache.roots[k][j];
					res = montmul64(pinv, 
							modsub64(qroot, proot,
							q2), q2, q2_w);

					if (res < lattice_size) {
						found_t *f = found_array + 
								my_threadid;
						f->p = p;
						f->q = q;
						f->which_poly = k;
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
