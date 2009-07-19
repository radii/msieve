/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 4/3/09
--------------------------------------------------------------------*/

#ifndef _GNFS_SQRT_SQRT_H_
#define _GNFS_SQRT_SQRT_H_

#include "gnfs.h"
#include <ap.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32 get_prime_for_sqrt(mp_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out); 

void alg_square_root(msieve_obj *obj, mp_poly_t *monic_alg_poly, 
			mp_t *n, mp_t *c, signed_mp_t *m1, signed_mp_t *m0,
			abpair_t *rlist, uint32 num_relations, uint32 check_q,
			mp_t *sqrt_a);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_SQRT_SQRT_H_ */
