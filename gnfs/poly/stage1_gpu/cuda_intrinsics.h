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

#if defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H)
#define CUDA_INTRINSICS_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef int int32;
typedef unsigned int uint32;
typedef unsigned long long uint64;

#define MIN(x, y) ((x) < (y) ? (x) : (y))

#define MASK24 0xffffff

/* 96-bit integers */

typedef struct {
	uint32 w[3];
} uint96;

/* 128-bit integers */

typedef struct {
	uint32 w[4];
} uint128;

/*------------------- Low-level functions ------------------------------*/

__device__ uint32
__uaddo(uint32 a, uint32 b) {
	uint32 res;
	asm("add.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ uint32
__uaddc(uint32 a, uint32 b) {
	uint32 res;
	asm("addc.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ uint32
__usubo(uint32 a, uint32 b) {
	uint32 res;
	asm("sub.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ uint32
__usubc(uint32 a, uint32 b) {
	uint32 res;
	asm("subc.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ uint32
__umul24hi(uint32 a, uint32 b) {
	uint32 res;
	asm("mul24.hi.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

/*------------------------ Comparison ---------------------------------*/

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

__device__ int32
cmp128(uint128 a, uint128 b)
{
	if (a.w[3] > b.w[3])
		return 1;
	if (a.w[3] < b.w[3])
		return -1;

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

/*-------------------- Addition/Subtraction --------------------------*/
__device__ uint96
add96(uint96 a, uint96 b)
{
	uint96 res;

	res.w[0] = __uaddo(a.w[0], b.w[0]);
	res.w[1] = __uaddc(a.w[1], b.w[1]);
	res.w[2] = __uaddc(a.w[2], b.w[2]);
	return res;
}

__device__ uint96
sub96(uint96 a, uint96 b)
{
	uint96 res;

	res.w[0] = __usubo(a.w[0], b.w[0]);
	res.w[1] = __usubc(a.w[1], b.w[1]);
	res.w[2] = __usubc(a.w[2], b.w[2]);
	return res;
}

__device__ uint128
add128(uint128 a, uint128 b)
{
	uint128 res;

	res.w[0] = __uaddo(a.w[0], b.w[0]);
	res.w[1] = __uaddc(a.w[1], b.w[1]);
	res.w[2] = __uaddc(a.w[2], b.w[2]);
	res.w[3] = __uaddc(a.w[3], b.w[3]);
	return res;
}

__device__ uint128
sub128(uint128 a, uint128 b)
{
	uint128 res;

	res.w[0] = __usubo(a.w[0], b.w[0]);
	res.w[1] = __usubc(a.w[1], b.w[1]);
	res.w[2] = __usubc(a.w[2], b.w[2]);
	res.w[3] = __usubc(a.w[3], b.w[3]);
	return res;
}

/*----------------- Squaring ----------------------------------------*/

__device__ uint64 
wide_sqr32(uint32 a)
{
	uint32 a0, a1;

	asm("{ .reg .u64 %dprod; \n\t"
	    "mul.wide.u32 %dprod, %2, %2; \n\t"
	    "cvt.u32.u64 %0, %dprod;      \n\t"
	    "shr.u64 %dprod, %dprod, 32;  \n\t"
	    "cvt.u32.u64 %1, %dprod;      \n\t"
	    "}                   \n\t"
	    : "=r"(a0), "=r"(a1)
	    : "r"(a));

	return (uint64)a1 << 32 | a0;
}

__device__ uint96 
wide_sqr48(uint64 a)
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

__device__ uint128
wide_sqr64(uint64 a)
{
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint64 acc;
	uint32 prod_lo, prod_hi;
	uint128 res;

	prod_lo = a0 * a0;
	prod_hi = __umulhi(a0, a0);
	res.w[0] = prod_lo;
	acc = (uint64)prod_hi;

	prod_lo = a0 * a1;
	prod_hi = __umulhi(a0, a1);
	acc = acc + prod_lo + prod_lo;
	res.w[1] = (uint32)acc;
	acc = (acc >> 32) + prod_hi + prod_hi;

	prod_lo = a1 * a1;
	prod_hi = __umulhi(a1, a1);
	acc += (uint64)prod_hi << 32 | prod_lo;
	res.w[2] = (uint32)acc;
	res.w[3] = (uint32)(acc >> 32);

	return res;
}

/* -------------------- Modular subtraction ------------------------*/

__device__ uint64 
modsub64(uint64 a, uint64 b, uint64 p) 
{
	uint32 r0, r1;
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 p0 = (uint32)p;
	uint32 p1 = (uint32)(p >> 32);

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %2, %4;        \n\t"
	    "subc.cc.u32 %1, %3, %5;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.cc.u32 %0, %0, %6; \n\t"
	    "@%pborrow addc.u32 %1, %1, %7; \n\t"
	    "} \n\t"
	    : "=r"(r0), "=r"(r1)
	    : "r"(a0), "r"(a1), 
	      "r"(b0), "r"(b1), 
	      "r"(p0), "r"(p1));

	return (uint64)r1 << 32 | r0;
}

__device__ uint96
modsub96(uint96 a, uint96 b, uint96 p) 
{
	uint96 res;
	uint32 a0 = a.w[0];
	uint32 a1 = a.w[1];
	uint32 a2 = a.w[2];
	uint32 b0 = b.w[0];
	uint32 b1 = b.w[1];
	uint32 b2 = b.w[2];
	uint32 p0 = p.w[0];
	uint32 p1 = p.w[1];
	uint32 p2 = p.w[2];

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %3, %6;        \n\t"
	    "subc.cc.u32 %1, %4, %7;        \n\t"
	    "subc.cc.u32 %2, %5, %8;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.cc.u32 %0, %0, %9; \n\t"
	    "@%pborrow addc.cc.u32 %1, %1, %10; \n\t"
	    "@%pborrow addc.u32 %2, %2, %11; \n\t"
	    "} \n\t"
	    : "=r"(res.w[0]), "=r"(res.w[1]), "=r"(res.w[2])
	    : "r"(a0), "r"(a1), "r"(a2),
	      "r"(b0), "r"(b1), "r"(b2),
	      "r"(p0), "r"(p1), "r"(p2));

	return res;
}

__device__ uint128 
modsub128(uint128 a, uint128 b, uint128 p) 
{
	uint128 res;
	uint32 a0 = a.w[0];
	uint32 a1 = a.w[1];
	uint32 a2 = a.w[2];
	uint32 a3 = a.w[3];
	uint32 b0 = b.w[0];
	uint32 b1 = b.w[1];
	uint32 b2 = b.w[2];
	uint32 b3 = b.w[3];
	uint32 p0 = p.w[0];
	uint32 p1 = p.w[1];
	uint32 p2 = p.w[2];
	uint32 p3 = p.w[3];

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %4, %8;        \n\t"
	    "subc.cc.u32 %1, %5, %9;        \n\t"
	    "subc.cc.u32 %2, %6, %10;        \n\t"
	    "subc.cc.u32 %3, %7, %11;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.cc.u32 %0, %0, %12; \n\t"
	    "@%pborrow addc.cc.u32 %1, %1, %13; \n\t"
	    "@%pborrow addc.cc.u32 %2, %2, %14; \n\t"
	    "@%pborrow addc.u32 %3, %3, %15; \n\t"
	    "} \n\t"
	    : "=r"(res.w[0]), "=r"(res.w[1]), "=r"(res.w[2]), "=r"(res.w[3])
	    : "r"(a0), "r"(a1), "r"(a2), "r"(a3),
	      "r"(b0), "r"(b1), "r"(b2), "r"(b3),
	      "r"(p0), "r"(p1), "r"(p2), "r"(p3));

	return res;
}

/*-------------------------- Modular inverse -------------------------*/

__device__ uint32 
modinv32(uint32 a, uint32 p) {

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

__device__ uint64 
modinv64(uint64 a, uint64 p) {

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

/*------------------- Montgomery arithmetic --------------------------*/
__device__ uint64 
montmul48(uint64 a, uint64 b,
		uint64 n, uint32 w) {

	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 24);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 24);
	uint32 n0 = (uint32)n;
	uint32 n1 = (uint32)(n >> 24);
	uint32 acc0, acc1;
	uint32 q0, q1;
	uint32 prod_lo, prod_hi;
	uint64 r;

	acc0 = __umul24(a0, b0);
	acc1 = __umul24hi(a0, b0) >> 16;
	q0 = __umul24(acc0, w);
	prod_lo = __umul24(q0, n0);
	prod_hi = __umul24hi(q0, n0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a0, b1);
	prod_hi = __umul24hi(a0, b1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(a1, b0);
	prod_hi = __umul24hi(a1, b0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q0, n1);
	prod_hi = __umul24hi(q0, n1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	q1 = __umul24(acc0, w);
	prod_lo = __umul24(q1, n0);
	prod_hi = __umul24hi(q1, n0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a1, b1);
	prod_hi = __umul24hi(a1, b1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(q1, n1);
	prod_hi = __umul24hi(q1, n1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);

	r = (uint64)acc1 << 32 | acc0;
	if (r >= n)
		return r - n;
	else
		return r;
}

__device__ uint64 
montmul64(uint64 a, uint64 b,
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

__device__ uint96 
montmul72(uint96 a, uint96 b,
		uint96 n, uint32 w) {

	uint32 a0 = a.w[0];
	uint32 a1 = a.w[0] >> 24 | a.w[1] << 8;
	uint32 a2 = a.w[1] >> 16 | a.w[2] << 16;
	uint32 b0 = b.w[0];
	uint32 b1 = b.w[0] >> 24 | b.w[1] << 8;
	uint32 b2 = b.w[1] >> 16 | b.w[2] << 16;
	uint32 n0 = n.w[0];
	uint32 n1 = n.w[0] >> 24 | n.w[1] << 8;
	uint32 n2 = n.w[1] >> 16 | n.w[2] << 16;
	uint32 acc0, acc1;
	uint32 q0, q1, q2;
	uint32 prod_lo, prod_hi;
	uint96 r;

	acc0 = __umul24(a0, b0);
	acc1 = __umul24hi(a0, b0) >> 16;
	q0 = __umul24(acc0, w);
	prod_lo = __umul24(q0, n0);
	prod_hi = __umul24hi(q0, n0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a0, b1);
	prod_hi = __umul24hi(a0, b1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(a1, b0);
	prod_hi = __umul24hi(a1, b0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q0, n1);
	prod_hi = __umul24hi(q0, n1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	q1 = __umul24(acc0, w);
	prod_lo = __umul24(q1, n0);
	prod_hi = __umul24hi(q1, n0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a0, b2);
	prod_hi = __umul24hi(a0, b2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(a1, b1);
	prod_hi = __umul24hi(a1, b1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(a2, b0);
	prod_hi = __umul24hi(a2, b0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q0, n2);
	prod_hi = __umul24hi(q0, n2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q1, n1);
	prod_hi = __umul24hi(q1, n1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	q2 = __umul24(acc0, w);
	prod_lo = __umul24(q2, n0);
	prod_hi = __umul24hi(q2, n0) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a1, b2);
	prod_hi = __umul24hi(a1, b2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(a2, b1);
	prod_hi = __umul24hi(a2, b1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q1, n2);
	prod_hi = __umul24hi(q1, n2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	prod_lo = __umul24(q2, n1);
	prod_hi = __umul24hi(q2, n1) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	r.w[0] = acc0 & MASK24;
	acc0 = acc0 >> 24 | acc1 << 8;

	prod_lo = __umul24(a2, b2);
	prod_hi = __umul24hi(a2, b2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(0, prod_hi);
	prod_lo = __umul24(q2, n2);
	prod_hi = __umul24hi(q2, n2) >> 16;
	acc0 = __uaddo(acc0, prod_lo);
	acc1 = __uaddc(acc1, prod_hi);
	r.w[0] |= acc0 << 24;
	r.w[1] = acc0 >> 8 | acc1 << 24;
	r.w[2] = acc1 >> 8;

	if (cmp96(r, n) >= 0)
		return sub96(r, n);
	else
		return r;
}

__device__ uint96 
montmul96(uint96 a, uint96 b,
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

__device__ uint128
montmul128(uint128 a, uint128 b,
		uint128 n, uint32 w) {

	uint32 acc0, acc1, acc2, acc3, acc4, nmult;
	uint32 prod_lo, prod_hi;
	uint64 prod;
	uint128 res;

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
	prod >>= 32;

	prod_lo = a.w[3] * b.w[0];
	prod_hi = __umulhi(a.w[3], b.w[0]);
	prod += (uint64)prod_hi << 32 | prod_lo;
	acc3 = (uint32)prod;
	acc4 = (uint32)(prod >> 32);

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

	prod_lo = nmult * n.w[3];
	prod_hi = __umulhi(nmult, n.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc2 = (uint32)prod;
	prod >>= 32;

	prod += acc4;
	acc3 = (uint32)prod;
	acc4 = (uint32)(prod >> 32);

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
	prod >>= 32;

	prod_lo = a.w[3] * b.w[1];
	prod_hi = __umulhi(a.w[3], b.w[1]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc3 = (uint32)prod;
	acc4 += (uint32)(prod >> 32);

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

	prod_lo = nmult * n.w[3];
	prod_hi = __umulhi(nmult, n.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc2 = (uint32)prod;
	prod >>= 32;

	prod += acc4;
	acc3 = (uint32)prod;
	acc4 = (uint32)(prod >> 32);

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
	prod >>= 32;

	prod_lo = a.w[3] * b.w[2];
	prod_hi = __umulhi(a.w[3], b.w[2]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc3 = (uint32)prod;
	acc4 += (uint32)(prod >> 32);

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

	prod_lo = nmult * n.w[3];
	prod_hi = __umulhi(nmult, n.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc2 = (uint32)prod;
	prod >>= 32;

	prod += acc4;
	acc3 = (uint32)prod;
	acc4 = (uint32)(prod >> 32);

	prod_lo = a.w[0] * b.w[3];   /*---------------------*/
	prod_hi = __umulhi(a.w[0], b.w[3]);
	prod = ((uint64)prod_hi << 32 | prod_lo) + acc0;
	acc0 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[1] * b.w[3];
	prod_hi = __umulhi(a.w[1], b.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc1;
	acc1 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[2] * b.w[3];
	prod_hi = __umulhi(a.w[2], b.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc2;
	acc2 = (uint32)prod;
	prod >>= 32;

	prod_lo = a.w[3] * b.w[3];
	prod_hi = __umulhi(a.w[3], b.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc3 = (uint32)prod;
	acc4 += (uint32)(prod >> 32);

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

	prod_lo = nmult * n.w[3];
	prod_hi = __umulhi(nmult, n.w[3]);
	prod += ((uint64)prod_hi << 32 | prod_lo) + acc3;
	acc2 = (uint32)prod;
	prod >>= 32;

	prod += acc4;
	acc3 = (uint32)prod;
	acc4 = (uint32)(prod >> 32);

	res.w[0] = acc0;        /*------------------------*/
	res.w[1] = acc1;
	res.w[2] = acc2;
	res.w[3] = acc3;
	if (acc4 > 0 || cmp128(res, n) >= 0)
		return sub128(res, n);
	else
		return res;
}

/*------------------ Initializing Montgomery arithmetic -----------------*/

__device__ uint32 
montmul24_w(uint32 n) {

	uint32 res = 8 - (n % 8);
	res = __umul24(res, 2 + __umul24(n, res));
	res = __umul24(res, 2 + __umul24(n, res));
	return __umul24(res, 2 + __umul24(n, res));
}

__device__ uint32 
montmul32_w(uint32 n) {

	uint32 res = 2 + n;
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	return res * (2 + n * res);
}

__device__ uint64 
montmul48_r(uint64 n, uint32 w) {

	uint32 shift;
	uint32 i;
	uint64 shifted_n;
	uint64 res;

	shift = __clzll(n);
	shifted_n = n << shift;
	res = -shifted_n;

	for (i = 64 - shift; i < 60; i++) {
		if (res >> 63)
			res = res + res - shifted_n;
		else
			res = res + res;

		if (res >= shifted_n)
			res -= shifted_n;
	}

	res = res >> shift;
	res = montmul48(res, res, n, w);
	return montmul48(res, res, n, w);
}

__device__ uint64 
montmul64_r(uint64 n, uint32 w) {

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
	res = montmul64(res, res, n, w);
	res = montmul64(res, res, n, w);
	return montmul64(res, res, n, w);
}

__device__ uint96 
montmul72_r(uint96 n, uint32 w) {

	/* 2^32 <= n < 2^72 */

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
	for (i = 95 - (word_shift + shift); i < 81; i++) {
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

	res = montmul72(res, res, n, w);
	res = montmul72(res, res, n, w);
	return montmul72(res, res, n, w);
}

__device__ uint96 
montmul96_r(uint96 n, uint32 w) {

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

	res = montmul96(res, res, n, w);
	res = montmul96(res, res, n, w);
	res = montmul96(res, res, n, w);
	return montmul96(res, res, n, w);
}

__device__ uint128
montmul128_r(uint128 n, uint32 w) {

	/* 2^64 <= n < 2^128 */

	uint32 shift, word_shift, comp_shift;
	uint32 i;
	uint128 shifted_n;
	uint128 res;

	if (n.w[3] == 0) {
		shifted_n.w[3] = n.w[2];
		shifted_n.w[2] = n.w[1];
		shifted_n.w[1] = n.w[0];
		shifted_n.w[0] = 0;
		word_shift = 32;
	}
	else {
		shifted_n = n;
		word_shift = 0;
	}

	shift = __clz(shifted_n.w[3]);
	comp_shift = 32 - shift;

	if (shift > 0) {
		shifted_n.w[3] = shifted_n.w[3] << shift | 
				shifted_n.w[2] >> comp_shift;
		shifted_n.w[2] = shifted_n.w[2] << shift | 
				shifted_n.w[1] >> comp_shift;
		shifted_n.w[1] = shifted_n.w[1] << shift | 
				shifted_n.w[0] >> comp_shift;
		shifted_n.w[0] = shifted_n.w[0] << shift;
	}

	res.w[0] = 0;
	res.w[1] = 0;
	res.w[2] = 0;
	res.w[3] = 0x80000000;
	for (i = 127 - (word_shift + shift); i < 136; i++) {
		if (res.w[3] & 0x80000000) {
			res = add128(res, res);
			res = sub128(res, shifted_n);
		}
		else {
			res = add128(res, res);
		}

		if (cmp128(res, shifted_n) > 0)
			res = sub128(res, shifted_n);
	}

	if (shift > 0) {
		res.w[0] = res.w[0] >> shift | res.w[1] << comp_shift;
		res.w[1] = res.w[1] >> shift | res.w[2] << comp_shift;
		res.w[2] = res.w[2] >> shift | res.w[3] << comp_shift;
		res.w[3] = res.w[3] >> shift;
	}
	if (word_shift > 0) {
		res.w[0] = res.w[1];
		res.w[1] = res.w[2];
		res.w[2] = res.w[3];
		res.w[3] = 0;
	}

	res = montmul128(res, res, n, w);
	res = montmul128(res, res, n, w);
	res = montmul128(res, res, n, w);
	return montmul128(res, res, n, w);
}

#ifdef __cplusplus
}
#endif

#endif /* defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H) */

