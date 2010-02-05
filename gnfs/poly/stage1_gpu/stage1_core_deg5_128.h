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

#ifndef _STAGE1_CORE_DEG5_128_H_
#define _STAGE1_CORE_DEG5_128_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __CUDACC__
	typedef int int32;
	typedef unsigned int uint32;
	typedef unsigned long long uint64;

	#define POLY_BATCH_SIZE 40

	/* 128-bit integers */

	typedef struct {
		uint32 w[4];
	} uint128;
#endif


/* structure indicating a collision */

typedef struct {
	uint64 p;
	uint64 q;
	uint32 which_poly;
	uint128 offset;
	uint128 proot;
} found_t;

#define P_SOA_BATCH_SIZE 2048

typedef struct {
	uint64 p[P_SOA_BATCH_SIZE];
	uint64 lattice_size[P_SOA_BATCH_SIZE];
	uint32 roots[4 * POLY_BATCH_SIZE][P_SOA_BATCH_SIZE];
} p_soa_t;

#define Q_SOA_BATCH_SIZE (5*30*256)

typedef struct {
	uint64 p[Q_SOA_BATCH_SIZE];
	uint32 roots[4 * (POLY_BATCH_SIZE + 1)][Q_SOA_BATCH_SIZE];
} q_soa_t;


#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_CORE_DEG5_128_H_ */
