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
	uint32 acc0, acc1, acc2, nmult;
	uint32 prod_lo, prod_hi;
	uint64 prod;

	prod_lo = a0 * b0;
	prod_hi = __umulhi(a0, b0);
	acc0 = prod_lo;

	prod = (uint64)prod_hi;
	prod_lo = a1 * b0;
	prod_hi = __umulhi(a1, b0);
	prod += ((uint64)prod_hi << 32 | prod_lo);
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	nmult = acc0 * w;

	prod_lo = nmult * n0;
	prod_hi = __umulhi(nmult, n0);
	prod = acc0 + ((uint64)prod_hi << 32 | prod_lo);
	prod = prod >> 32;

	prod_lo = nmult * n1;
	prod_hi = __umulhi(nmult, n1);
	prod += (uint64)acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = (prod >> 32) + (uint64)acc2;
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	prod_lo = a0 * b1;
	prod_hi = __umulhi(a0, b1);
	prod = (uint64)acc0 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = prod >> 32;

	prod_lo = a1 * b1;
	prod_hi = __umulhi(a1, b1);
	prod += (uint64)acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32) + acc2;

	nmult = acc0 * w;

	prod_hi = __umulhi(nmult, n0);
	prod_lo = nmult * n0;
	prod = acc0 + ((uint64)prod_hi << 32 | prod_lo);
	prod = prod >> 32;

	prod_hi = __umulhi(nmult, n1);
	prod_lo = nmult * n1;
	prod += acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = (prod >> 32) + (uint64)acc2;
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	prod = (uint64)acc1 << 32 | acc0;
	if (acc2 || prod >= n)
		return prod - n;
	else
		return prod;
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
__device__ p_packed_t *
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + curr->num_roots);
}

/*------------------------------------------------------------------------*/
__constant__ uint64 pbatch[P_ARRAY_WORDS];

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
		uint64 q2 = (uint64)q * q;
		uint32 q2_w = montmul_w((uint32)q2);
		uint64 q2_r = montmul_r(q2, q2_w);
		p_packed_t *curr_p = (p_packed_t *)pbatch;
		
		for (j = 0; j < num_p; j++) {
			uint32 p = curr_p->p;
			uint64 p2 = (uint64)p * p;
			uint32 pinvmodq = modinv(p, q);

			uint32 num_proots = curr_p->num_roots;
			uint32 lattice_size = curr_p->lattice_size;
			uint64 pinv, tmp;

			tmp = (uint64)pinvmodq * pinvmodq;
			tmp = montmul(tmp, q2_r, q2, q2_w);
			pinv = montmul(p2, tmp, q2, q2_w);
			pinv = modsub((uint64)2, pinv, q2);
			pinv = montmul(pinv, tmp, q2, q2_w);
			pinv = montmul(pinv, q2_r, q2, q2_w);

			for (k = 0; k < num_qroots; k++) {

				uint64 qroot = qbatch->roots[k][i];

				for (m = 0; m < num_proots; m++) {

					uint64 proot = curr_p->roots[m];
					uint64 res = montmul(pinv, 
							modsub(qroot, 
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
