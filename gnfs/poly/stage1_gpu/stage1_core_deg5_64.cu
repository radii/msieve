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

#include "cuda_intrinsics.h"
#include "stage1_core_deg5_64.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MIN(x, y) ((x) < (y) ? (x) : (y))

/*------------------------------------------------------------------------*/
__device__ uint64 
modsub(uint64 a, uint64 b, uint64 p) 
{
	uint64 t = 0, tr;
	tr = a - b;
	if (tr > a)
		t = p;
	return tr + t;
}

/*------------------------------------------------------------------------*/
__device__ uint32 
modinv(uint32 a, uint32 p) {

	uint32 ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend - q * divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

/*------------------------------------------------------------------------*/
__device__ uint32 
montmul_w(uint32 n) {

	uint32 res = 2 + n;
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	return res * (2 + n * res);
}

/*------------------------------------------------------------------------*/
__device__ uint64 
montmul(uint64 a, uint64 b,
		uint64 n, uint32 w) {

	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 n0 = (uint32)n;
	uint32 n1 = (uint32)(n >> 32);
	uint32 acc0, acc1, acc2;
	uint32 q0, q1;
	uint32 prod_lo, prod_hi;
	uint64 r;

	acc0 = a0 * b0;
	acc1 = __umulhi(a0, b0);
	q0 = acc0 * w;
	prod_lo = q0 * n0;
	prod_hi = __umulhi(q0, n0);
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc2 = __uaddc(0, 0);

	prod_lo = a0 * b1;
	prod_hi = __umulhi(a0, b1);
	acc0 = __uaddo(acc1, prod_lo);
	acc1 = __uaddc(acc2, prod_hi);
	acc2 = __uaddc(0, 0);
	prod_lo = a1 * b0;
	prod_hi = __umulhi(a1, b0);
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc2 = __uaddc(acc2, 0);
	prod_lo = q0 * n1;
	prod_hi = __umulhi(q0, n1);
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc2 = __uaddc(acc2, 0);
	q1 = acc0 * w;
	prod_lo = q1 * n0;
	prod_hi = __umulhi(q1, n0);
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc2 = __uaddc(acc2, 0);

	prod_lo = a1 * b1;
	prod_hi = __umulhi(a1, b1);
	acc0 = __uaddo(acc1, prod_lo);
	acc1 = __uaddc(acc2, prod_hi);
	acc2 = __uaddc(0, 0);
	prod_lo = q1 * n1;
	prod_hi = __umulhi(q1, n1);
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc2 = __uaddc(acc2, 0);

	r = (uint64)acc1 << 32 | acc0;
	if (acc2 || r >= n)
		return r - n;
	else
		return r;
}

/*------------------------------------------------------------------------*/
__device__ uint64 
montmul_r(uint64 n, uint32 w) {

	uint32 shift;
	uint32 i;
	uint64 shifted_n;
	uint64 res;

	shift = __clzll(n);
	shifted_n = n << shift;
	res = -shifted_n;

	for (i = 64 - shift; i < 72; i++) {
		if (res >> 63)
			res = res + res - shifted_n;
		else
			res = res + res;

		if (res >= shifted_n)
			res -= shifted_n;
	}

	res = res >> shift;
	res = montmul(res, res, n, w);
	res = montmul(res, res, n, w);
	return montmul(res, res, n, w);
}

/*------------------------------------------------------------------------*/
#define SHARED_BATCH_SIZE 48

typedef struct {
	uint32 p[SHARED_BATCH_SIZE];
	uint32 lattice_size[SHARED_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE][SHARED_BATCH_SIZE];
} p_soa_shared_t;

__shared__ p_soa_shared_t pbatch_cache;

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
		q2_w = montmul_w((uint32)q2);
		q2_r = montmul_r(q2, q2_w);

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
				uint32 pinvmodq = modinv(p, q);

				uint32 lattice_size = 
						pbatch_cache.lattice_size[j];
				uint64 pinv, tmp;

				tmp = (uint64)pinvmodq * pinvmodq;
				tmp = montmul(tmp, q2_r, q2, q2_w);
				pinv = montmul(p2, tmp, q2, q2_w);
				pinv = modsub((uint64)2, pinv, q2);
				pinv = montmul(pinv, tmp, q2, q2_w);
				pinv = montmul(pinv, q2_r, q2, q2_w);

				for (k = 0; k < num_roots; k++) {

					uint64 qroot;
					uint64 proot;
					uint64 res;

					qroot = prefetch;
					prefetch = qbatch->roots[k+1][i];
					proot = pbatch_cache.roots[k][j];
					res = montmul(pinv, modsub(qroot, proot,
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
