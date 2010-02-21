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

#ifndef _STAGE1_CORE_DEG5_64_H_
#define _STAGE1_CORE_DEG5_64_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __CUDACC__
	#include "cuda_intrinsics.h"

	#define POLY_BATCH_SIZE 40
#endif

/* structure indicating a collision */

typedef struct {
	uint32 p;
	uint32 q;
	uint32 which_poly;
	uint32 pad;
	uint64 offset;
	uint64 proot;
} found_t;

/* the outer loop needs parallel access to different p,
   so we store in SOA format. All the entries in the structure
   have the same number of roots */

#define P_SOA_BATCH_SIZE 2048

typedef struct {
	uint32 p[P_SOA_BATCH_SIZE];
	uint32 lattice_size[P_SOA_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE][P_SOA_BATCH_SIZE];
} p_soa_t;

#define Q_SOA_BATCH_SIZE (3*30*384)

typedef struct {
	uint32 p[Q_SOA_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE+1][Q_SOA_BATCH_SIZE];
} q_soa_t;


#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_CORE_DEG5_64_H_ */
