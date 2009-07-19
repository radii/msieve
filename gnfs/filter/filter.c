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

#include "filter.h"

/*--------------------------------------------------------------------*/
static void find_fb_size(factor_base_t *fb, 
			uint32 limit_r, uint32 limit_a,
			uint32 *entries_r_out, uint32 *entries_a_out) {

	prime_sieve_t prime_sieve;
	uint32 entries_r = 0;
	uint32 entries_a = 0;

	init_prime_sieve(&prime_sieve, 0, 
			MAX(limit_r, limit_a) + 1000);

	while (1) {
		uint32 p = get_next_prime(&prime_sieve);
		uint32 num_roots;
		uint32 high_coeff;
		uint32 roots[MAX_POLY_DEGREE + 1];

		if (p >= limit_r && p >= limit_a)
			break;

		if (p < limit_r) {
			num_roots = poly_get_zeros(roots, &fb->rfb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_r += num_roots;
		}

		if (p < limit_a) {
			num_roots = poly_get_zeros(roots, &fb->afb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_a += num_roots;
		}
	}

	free_prime_sieve(&prime_sieve);
	*entries_r_out = entries_r;
	*entries_a_out = entries_a;
}

/*--------------------------------------------------------------------*/
#define MAX_KEEP_WEIGHT 45

/* the multiple of the amount of excess needed for
   merging to proceed */

#define FINAL_EXCESS_FRACTION 1.16

uint32 nfs_filter_relations(msieve_obj *obj, mp_t *n) {

	filter_t filter;
	merge_t merge;
	uint32 filtmin_r;
	uint32 filtmin_a;
	uint32 entries_r;
	uint32 entries_a;
	uint32 extra_needed;
	uint32 max_weight;
	factor_base_t fb;
	double cpu_time;

	logprintf(obj, "\n");
	logprintf(obj, "commencing relation filtering\n");

	cpu_time = get_cpu_time();
	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly, NULL)) {
		printf("filtering failed to read polynomials\n");
		exit(-1);
	}

	/* delete duplicate relations, and determine the cutoff
	   size of large primes used in the rest of the filtering.
	   We do not just use the factor base, since it may be too 
	   small or too large for this particular dataset */

	filtmin_r = filtmin_a = nfs_purge_duplicates(obj, &fb, 
						(uint32)obj->nfs_upper);
	if (obj->nfs_lower)
		filtmin_r = filtmin_a = obj->nfs_lower;

	/* singleton removal happens in at least two phases. 
	   Phase 1 does the singleton removal for all the 
	   ideals larger than filtmin, regardless of ideal weight.
	   
	   Phase 2 always runs, even if there is not enough excess 
	   for phase 1; the first pass is just to remove the 
	   majority of the singletons */

	memset(&filter, 0, sizeof(filter));
	filter.filtmin_r = filtmin_r;
	filter.filtmin_a = filtmin_a;
	logprintf(obj, "reading rational ideals above %u\n", filtmin_r);
	logprintf(obj, "reading algebraic ideals above %u\n", filtmin_a);

	nfs_purge_singletons_initial(obj, &fb, &filter);

	/* throw away the current relation list */

	free(filter.relation_array);
	filter.relation_array = NULL;

	/* perform the other filtering passes; these use a very small
	   filtering bound, and rely on most unneeded relations
	   having been deleted already. The set of relations remaining
	   will be forwarded to the final merge phase */

	if (filtmin_r < 3000000 && filtmin_a < 3000000) {
		/* small problems are filtered aggressively */
		filtmin_r = MIN(filtmin_r / 2, 100000);
		filtmin_a = MIN(filtmin_a / 2, 100000);
		max_weight = 25;
	}
	else {
		filtmin_r = MIN(filtmin_r / 2, 720000);
		filtmin_a = MIN(filtmin_a / 2, 720000);
		max_weight = 20;
	}

	for (; max_weight < MAX_KEEP_WEIGHT; max_weight += 5) {

		find_fb_size(&fb, filtmin_r, filtmin_a, &entries_r, &entries_a);
		filter.filtmin_r = filtmin_r;
		filter.filtmin_a = filtmin_a;
		filter.target_excess = entries_r + entries_a;
		logprintf(obj, "reading rational ideals above %u\n", 
						filter.filtmin_r);
		logprintf(obj, "reading algebraic ideals above %u\n", 
						filter.filtmin_a);

		nfs_purge_singletons(obj, &fb, &filter, max_weight);

		/* make the clique removal more conservative by
		   leaving some of the excess; this makes the merge 
		   phase easier. Note that the singleton removal
		   probably threw away many large ideals that occur
		   too often to be worth tracking, which forces the 
		   target matrix size to increase, so that target_excess 
		   is larger now */

		extra_needed = filter.target_excess;
		filter.target_excess *= FINAL_EXCESS_FRACTION;

		/* give up if there is not enough excess 
		   to form a matrix */

		if (filter.num_relations < filter.num_ideals ||
		    filter.num_relations - filter.num_ideals < 
		    				filter.target_excess) {
			uint32 relations_needed = 1000000;

			if (filter.num_relations > filter.num_ideals) {
				relations_needed = 3 * (filter.target_excess -
						(filter.num_relations - 
					 	 filter.num_ideals));
				relations_needed = MAX(relations_needed, 20000);
			}
			free(filter.relation_array);
			filter.relation_array = NULL;
			return relations_needed;
		}

		merge.num_extra_relations = NUM_EXTRA_RELATIONS;

		/* build the matrix */

		filter_make_relsets(obj, &filter, &merge, extra_needed);

		/* retry the filtering with more aggressive
		   settings if the clique max weight was small and
		   the largest number of relations in a relation set
		   is close to the clique max weight. This is an
		   indication that the clique processing ignored too
		   many ideals that the merge phase should have seen.
		   We do not preemptively rerun the clique processing 
		   because for large problems rereading the dataset is 
		   much more expensive than just redoing the merge */

		if (filter.max_ideal_weight < 18 &&
		    merge.max_relations > filter.max_ideal_weight - 3) {
			logprintf(obj, "matrix can improve, retrying\n");
			filter_free_relsets(&merge);
			continue;
		}

		/* do not accept the collection of generated cycles
		   unless the matrix they form is dense enough. It's
		   unfortunate, but different datasets could look the
		   same going into the merge phase but lead to very
		   different matrix densities. The problem is that if
		   the matrix is not dense enough experience shows
		   that the linear algebra sometimes is unable to find
		   nontrivial dependencies. This is especially common
		   when the bound on large primes is much larger than
		   is sensible. In fact, I haven't seen a matrix
		   with average cycle weight under 50.0 get solved
		   successfully */

		if (merge.avg_cycle_weight > 63.0)
			break;

		logprintf(obj, "matrix not dense enough, retrying\n");
		filter_free_relsets(&merge);

		/* lower the filtering bound for next time */

		filtmin_r = 0.8 * filtmin_r;
		filtmin_a = 0.8 * filtmin_a;
	}

	if (max_weight >= MAX_KEEP_WEIGHT) {
		printf("error: too many merge attempts\n");
		exit(-1);
	}

	/* optimize the collection of cycles 
	   and save the result */

	filter_postproc_relsets(obj, &merge);
	filter_dump_relsets(obj, &merge);
	filter_free_relsets(&merge);
	cpu_time = get_cpu_time() - cpu_time;
	logprintf(obj, "RelProcTime: %u\n", (uint32)cpu_time);
	return 0;
}
