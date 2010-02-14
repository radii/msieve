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

#include "stage1_core_deg46_64.h"

#ifdef __cplusplus
extern "C" {
#endif

__constant__ uint64 pbatch[P_ARRAY_WORDS];

/*------------------------------------------------------------------------*/
__device__ p_packed_t *
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + curr->num_roots);
}

/*------------------------------------------------------------------------*/
__global__ void
sieve_kernel_64(q_soa_t *qbatch, 
             uint32 num_q,
	     uint32 num_qroots,
	     uint32 num_p,
	     found_t *found_array)
{
	uint32 my_threadid;
	uint32 num_threads;
	uint32 i, j, k, m;

	my_threadid = blockIdx.x * blockDim.x + threadIdx.x;
	num_threads = gridDim.x * blockDim.x;
	found_array[my_threadid].p = 0;

	for (i = my_threadid; i < num_q; i += num_threads) {
		uint32 q = qbatch->p[i];
		uint64 q2 = wide_sqr32(q);
		uint32 q2_w = montmul32_w((uint32)q2);
		uint64 q2_r = montmul64_r(q2, q2_w);
		p_packed_t *curr_p = (p_packed_t *)pbatch;
		
		for (j = 0; j < num_p; j++) {
			uint32 p = curr_p->p;
			uint64 p2 = wide_sqr32(p);
			uint32 pinvmodq = modinv32(p, q);

			uint32 num_proots = curr_p->num_roots;
			uint32 lattice_size = curr_p->lattice_size;
			uint64 pinv, tmp;

			tmp = wide_sqr32(pinvmodq);
			tmp = montmul64(tmp, q2_r, q2, q2_w);
			pinv = montmul64(p2, tmp, q2, q2_w);
			pinv = modsub64((uint64)2, pinv, q2);
			pinv = montmul64(pinv, tmp, q2, q2_w);
			pinv = montmul64(pinv, q2_r, q2, q2_w);

			for (k = 0; k < num_qroots; k++) {

				uint64 qroot = qbatch->roots[k][i];

				for (m = 0; m < num_proots; m++) {

					uint64 proot = curr_p->roots[m];
					uint64 res = montmul64(pinv, 
							modsub64(qroot, 
								proot, q2),
							q2, q2_w);

					if (res < lattice_size) {
						found_t *f = found_array +
								my_threadid;
						f->p = p;
						f->q = q;
						f->offset = res;
						f->proot = proot;
					}
				}
			}

			curr_p = p_packed_next(curr_p);
		}
	}
}

#ifdef __cplusplus
}
#endif
