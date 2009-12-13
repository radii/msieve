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

/*------------------------------------------------------------------------*/
__device__ int32
cmp96(uint96 a, uint96 b)
{
	if (a.w[2] > b.w[2])
		return 1;
	if (a.w[2] < b.w[2])
		return -1;

	if (a.w[1] > b.w[1])
		return 1;
	if (a.w[1] < b.w[1])
		return -1;

	if (a.w[0] > b.w[0])
		return 1;
	if (a.w[0] < b.w[0])
		return -1;
	return 0;
}

/*------------------------------------------------------------------------*/
__device__ uint96
add96(uint96 a, uint96 b)
{
	uint32 c;
	uint32 acc;
	uint96 res;

	acc = a.w[0] + b.w[0];
	res.w[0] = acc;
	c = (acc < a.w[0]);

	acc = a.w[1] + c;
	c = (acc < a.w[1]);
	res.w[1] = acc + b.w[1];
	c += (res.w[1] < acc);

	res.w[2] = a.w[2] + b.w[2] + c;
	return res;
}

/*------------------------------------------------------------------------*/
__device__ uint96
sub96(uint96 a, uint96 b)
{
	uint32 c;
	uint32 acc;
	uint96 res;

	acc = a.w[0] - b.w[0];
	res.w[0] = acc;
	c = (acc > a.w[0]);

	acc = a.w[1] - c;
	c = (acc > a.w[1]);
	res.w[1] = acc - b.w[1];
	c += (res.w[1] > acc);

	res.w[2] = a.w[2] - b.w[2] - c;
	return res;
}

/*------------------------------------------------------------------------*/
__device__ uint96 
modsub(uint96 a, uint96 b, uint96 p) 
{
	/* this could be 7 branch-less instructions
	   if nvcc allowed inline asm */

	uint96 res = sub96(a, b);

	if (cmp96(res, a) > 0)
		res = add96(res, p);

	return res;
}

/*------------------------------------------------------------------------*/
__device__ uint96 
wide_sqr(uint64 a)
{
	/* a < 2^48 */

	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint64 acc;
	uint32 prod_lo, prod_hi;
	uint96 res;

	prod_lo = a0 * a0;
	prod_hi = __umulhi(a0, a0);
	res.w[0] = prod_lo;
	acc = (uint64)prod_hi;

	prod_lo = a0 * a1;
	prod_hi = __umulhi(a0, a1);
	acc += 2 * ((uint64)prod_hi << 32 | prod_lo);
	res.w[1] = (uint32)acc;
	res.w[2] = (uint32)(acc >> 32) + __umul24(a1, a1);

	return res;
}

/*------------------------------------------------------------------------*/
__device__ uint64 
modinv(uint64 a, uint64 p) {

	uint64 ps1, ps2, dividend, divisor, rem, q, t;
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
__device__ uint96 
montmul(uint96 a, uint96 b,
		uint96 n, uint32 w) {

	uint32 acc0, acc1, acc2, acc3, nmult;
	uint32 prod_lo, prod_hi;
	uint64 prod;
	uint96 res;

	acc0 = a.w[0] * b.w[0];   /*---------------------*/
	prod = (uint64)(__umulhi(a.w[0], b.w[0]));

	prod_lo = a.w[1] * b.w[0];
	prod_hi = __umulhi(a.w[1], b.w[0]);
	prod += (uint64)prod_hi << 32 | prod_lo;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[2] * b.w[0];
	prod_hi = __umulhi(a.w[2], b.w[0]);
	prod += (uint64)prod_hi << 32 | prod_lo;
	acc2 = (uint32)prod;
	acc3 = (uint32)(prod >> 32);

	nmult = acc0 * w;      /*------------------------*/

	prod_lo = nmult * n.w[0];
	prod_hi = __umulhi(nmult, n.w[0]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	prod >>= 32;

	prod_lo = nmult * n.w[1];
	prod_hi = __umulhi(nmult, n.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = nmult * n.w[2];
	prod_hi = __umulhi(nmult, n.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod += acc3;
	acc2 = (uint32)prod;
	acc3 = (uint32)(prod >> 32);

	prod_lo = a.w[0] * b.w[1];   /*---------------------*/
	prod_hi = __umulhi(a.w[0], b.w[1]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[1] * b.w[1];
	prod_hi = __umulhi(a.w[1], b.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[2] * b.w[1];
	prod_hi = __umulhi(a.w[2], b.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc2 = (uint32)prod;
	acc3 += (uint32)(prod >> 32);

	nmult = acc0 * w;      /*------------------------*/

	prod_lo = nmult * n.w[0];
	prod_hi = __umulhi(nmult, n.w[0]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	prod >>= 32;

	prod_lo = nmult * n.w[1];
	prod_hi = __umulhi(nmult, n.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = nmult * n.w[2];
	prod_hi = __umulhi(nmult, n.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod += acc3;
	acc2 = (uint32)prod;
	acc3 = (uint32)(prod >> 32);

	prod_lo = a.w[0] * b.w[2];   /*---------------------*/
	prod_hi = __umulhi(a.w[0], b.w[2]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[1] * b.w[2];
	prod_hi = __umulhi(a.w[1], b.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[2] * b.w[2];
	prod_hi = __umulhi(a.w[2], b.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc2 = (uint32)prod;
	acc3 += (uint32)(prod >> 32);

	nmult = acc0 * w;      /*------------------------*/

	prod_lo = nmult * n.w[0];
	prod_hi = __umulhi(nmult, n.w[0]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	prod >>= 32;

	prod_lo = nmult * n.w[1];
	prod_hi = __umulhi(nmult, n.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = nmult * n.w[2];
	prod_hi = __umulhi(nmult, n.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod += acc3;
	acc2 = (uint32)prod;
	acc3 = (uint32)(prod >> 32);

	res.w[0] = acc0;        /*------------------------*/
	res.w[1] = acc1;
	res.w[2] = acc2;
	if (acc3 > 0 || cmp96(res, n) >= 0)
		return sub96(res, n);
	else
		return res;
}

/*------------------------------------------------------------------------*/
__device__ uint96 
montmul_r(uint96 n, uint32 w) {

	/* 2^32 <= n < 2^96 */

	uint32 shift, word_shift, comp_shift;
	uint32 i;
	uint96 shifted_n;
	uint96 res;

	if (n.w[2] == 0) {
		shifted_n.w[2] = n.w[1];
		shifted_n.w[1] = n.w[0];
		shifted_n.w[0] = 0;
		word_shift = 32;
	}
	else {
		shifted_n = n;
		word_shift = 0;
	}

	shift = __clz(shifted_n.w[2]);
	comp_shift = 32 - shift;

	if (shift > 0) {
		shifted_n.w[2] = shifted_n.w[2] << shift | 
				shifted_n.w[1] >> comp_shift;
		shifted_n.w[1] = shifted_n.w[1] << shift | 
				shifted_n.w[0] >> comp_shift;
		shifted_n.w[0] = shifted_n.w[0] << shift;
	}

	res.w[0] = 0;
	res.w[1] = 0;
	res.w[2] = 0x80000000;
	for (i = 95 - (word_shift + shift); i < 102; i++) {
		if (res.w[2] & 0x80000000) {
			res = add96(res, res);
			res = sub96(res, shifted_n);
		}
		else {
			res = add96(res, res);
		}

		if (cmp96(res, shifted_n) > 0)
			res = sub96(res, shifted_n);
	}

	if (shift > 0) {
		res.w[0] = res.w[0] >> shift | res.w[1] << comp_shift;
		res.w[1] = res.w[1] >> shift | res.w[2] << comp_shift;
		res.w[2] = res.w[2] >> shift;
	}
	if (word_shift > 0) {
		res.w[0] = res.w[1];
		res.w[1] = res.w[2];
		res.w[2] = 0;
	}

	res = montmul(res, res, n, w);
	res = montmul(res, res, n, w);
	res = montmul(res, res, n, w);
	return montmul(res, res, n, w);
}

/*------------------------------------------------------------------------*/
__device__ p_packed_t *
p_packed_next(p_packed_t *curr)
{
	return (p_packed_t *)((uint64 *)curr + 
			P_PACKED_HEADER_WORDS + 
			3 * (curr->num_roots / 2));
}

/*------------------------------------------------------------------------*/
__constant__ uint64 pbatch[P_ARRAY_WORDS];

__constant__ uint96 two = {{2, 0, 0}};

__global__ void
sieve_kernel_96(q_soa_t *qbatch, 
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
		uint96 q2 = wide_sqr(q);
		uint32 q2_w = montmul_w(q2.w[0]);
		uint96 q2_r = montmul_r(q2, q2_w);
		p_packed_t *curr_p = (p_packed_t *)pbatch;
		
		for (j = 0; j < num_p; j++) {
			uint64 p = curr_p->p;
			uint96 p2 = wide_sqr(p);
			uint64 pinvmodq = modinv(p, q);

			uint32 num_proots = curr_p->num_roots;
			uint64 lattice_size = curr_p->lattice_size;
			uint96 pinv, tmp;
			uint96 test1;

			test1.w[0] = (uint32)lattice_size;
			test1.w[1] = (uint32)(lattice_size >> 32);
			test1.w[2] = 0;

			tmp = wide_sqr(pinvmodq);
			tmp = montmul(tmp, q2_r, q2, q2_w);
			pinv = montmul(p2, tmp, q2, q2_w);
			pinv = modsub(two, pinv, q2);
			pinv = montmul(pinv, tmp, q2, q2_w);
			pinv = montmul(pinv, q2_r, q2, q2_w);

			for (k = 0; k < 3 * num_qroots; k += 3) {

				uint96 qroot;

				qroot.w[0] = qbatch->roots[k][i];
				qroot.w[1] = qbatch->roots[k+1][i];
				qroot.w[2] = qbatch->roots[k+2][i];

				for (m = 0; m < num_proots; m++) {

					uint96 proot = curr_p->roots[m];
					uint96 res = montmul(pinv, 
							modsub(qroot, 
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
