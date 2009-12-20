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

#ifndef _STAGE1_CORE_DEG46_64_H_
#define _STAGE1_CORE_DEG46_64_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __CUDACC__
typedef unsigned int uint32;
typedef unsigned long long uint64;
#define MAX_ROOTS 128
#endif

/* structure indicating a collision */

typedef struct {
	uint32 p;
	uint32 q;
	uint64 offset;
	uint64 proot;
} found_t;

#define P_ARRAY_WORDS 1000

#define P_PACKED_HEADER_WORDS 2

typedef struct {
	uint32 p;
	uint32 lattice_size;
	uint32 num_roots;
	uint32 pad;
	uint64 roots[MAX_ROOTS];
} p_packed_t;

#define Q_SOA_BATCH_SIZE (3*30*384)

typedef struct {
	uint32 p[Q_SOA_BATCH_SIZE];
	uint64 roots[MAX_ROOTS][Q_SOA_BATCH_SIZE];
} q_soa_t;

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_CORE_DEG46_64_H_ */
