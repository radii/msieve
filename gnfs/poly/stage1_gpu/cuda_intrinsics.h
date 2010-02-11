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

__device__ unsigned int
__uaddo(unsigned int a, unsigned int b) {
	unsigned int res;
	asm("add.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ unsigned int
__uaddc(unsigned int a, unsigned int b) {
	unsigned int res;
	asm("addc.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ unsigned int
__usubo(unsigned int a, unsigned int b) {
	unsigned int res;
	asm("sub.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ unsigned int
__usubc(unsigned int a, unsigned int b) {
	unsigned int res;
	asm("subc.cc.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}

__device__ unsigned int
__umul24hi(unsigned int a, unsigned int b) {
	unsigned int res;
	asm("mul24.hi.u32 %0, %1, %2; /* inline */ \n\t" 
	    : "=r" (res) : "r" (a) , "r" (b));
	return res;
}


#ifdef __cplusplus
}
#endif

#endif

