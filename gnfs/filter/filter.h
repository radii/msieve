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

/* implementation of number field sieve filtering */

#ifndef _GNFS_FILTER_FILTER_H_
#define _GNFS_FILTER_FILTER_H_

#include <common/filter/filter.h>
#include "gnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* create '<savefile_name>.d', a binary file containing
   the line numbers of unique relations. The return value is
   the large prime bound to use for the rest of the filtering.
   Duplicate removal only applies to the first max_relations
   relations found (or all relations if zero) */

uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint32 max_relations); 

/* read '<savefile_name>.d' and create '<savefile_name>.s', a 
   binary file containing the line numbers of relations that
   are not singletons. All ideals larger than the bounds specified
   in 'filter' are tacked */
   
void nfs_purge_singletons_initial(msieve_obj *obj, 
			factor_base_t *fb, filter_t *filter);

/* read and modify '<savefile_name>.s', to account for ideals
   larger than the filtering bound in 'filter' and occurring in
   at most max_ideal_weight relations (or all relations if zero) */
   
void nfs_purge_singletons(msieve_obj *obj, factor_base_t *fb,
			filter_t *filter, uint32 max_ideal_weight);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_FILTER_FILTER_H_ */
