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

#include "stage1_core_deg6_96.h"

#ifdef __cplusplus
extern "C" {
#endif

__constant__ uint64 pbatch[P_ARRAY_WORDS];

__constant__ uint96 two = {{2, 0, 0}};

/*------------------------------------------------------------------------*/
__device__ p_packed_t *
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + 
			3 * (curr->num_roots / 2));
}

/*------------------------------------------------------------------------*/
__global__ void
sieve_kernel(q_soa_t *qbatch, 
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
		uint64 q = qbatch->p[i];
		uint96 q2 = wide_sqr48(q);
		uint32 q2_w = montmul24_w(q2.w[0]);
		uint96 q2_r = montmul72_r(q2, q2_w);
		p_packed_t *curr_p = (p_packed_t *)pbatch;
		
		for (j = 0; j < num_p; j++) {
			uint64 p = curr_p->p;
			uint96 p2 = wide_sqr48(p);
			uint64 pinvmodq = modinv64(p, q);

			uint32 num_proots = curr_p->num_roots;
			uint64 lattice_size = curr_p->lattice_size;
			uint96 pinv, tmp;
			uint96 test1;

			test1.w[0] = (uint32)lattice_size;
			test1.w[1] = (uint32)(lattice_size >> 32);
			test1.w[2] = 0;

			tmp = wide_sqr48(pinvmodq);
			tmp = montmul72(tmp, q2_r, q2, q2_w);
			pinv = montmul72(p2, tmp, q2, q2_w);
			pinv = modsub96(two, pinv, q2);
			pinv = montmul72(pinv, tmp, q2, q2_w);
			pinv = montmul72(pinv, q2_r, q2, q2_w);

			for (k = 0; k < 3 * num_qroots; k += 3) {

				uint96 qroot;

				qroot.w[0] = qbatch->roots[k][i];
				qroot.w[1] = qbatch->roots[k+1][i];
				qroot.w[2] = qbatch->roots[k+2][i];

				for (m = 0; m < num_proots; m++) {

					uint96 proot = curr_p->roots[m];
					uint96 res = montmul72(pinv, 
							modsub96(qroot, 
								proot, q2),
							q2, q2_w);

					if (cmp96(res, test1) <= 0) {
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
