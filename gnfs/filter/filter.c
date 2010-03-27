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

#include "filter.h"

/*--------------------------------------------------------------------*/
static void find_fb_size(factor_base_t *fb, 
			uint32 limit_r, uint32 limit_a,
			uint32 *entries_r_out, uint32 *entries_a_out) {

	prime_sieve_t prime_sieve;
	uint32 entries_r = 0;
	uint32 entries_a = 0;

	/* If the filtering bounds are extremely large, just 
	   estimate the target as the number of primes less than 
	   the filtering bounds */

	if (limit_r > 20000000 && limit_a > 20000000) {
		*entries_r_out = 1.02 * limit_r / (log((double)limit_r) - 1);
		*entries_a_out = 1.02 * limit_a / (log((double)limit_a) - 1);
		return;
	}

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
static uint32 check_excess(filter_t *filter) {

	/* give up if there is not enough excess 
	   to form a matrix */

	uint32 relations_needed = 0;

	if (filter->num_relations < filter->num_ideals ||
	    filter->num_relations - filter->num_ideals < 
	    				filter->target_excess) {
		relations_needed = 1000000;

		if (filter->num_relations > filter->num_ideals) {
			relations_needed = 3 * (filter->target_excess -
					(filter->num_relations - 
				 	 filter->num_ideals));
			relations_needed = MAX(relations_needed, 1000000);
		}
		free(filter->relation_array);
		filter->relation_array = NULL;
	}
	return relations_needed;
}

/*--------------------------------------------------------------------*/
static void dump_relation_numbers(msieve_obj *obj, filter_t *filter) {

	uint32 i;
	char buf[256];
	FILE *relation_fp;
	relation_ideal_t *r = filter->relation_array;

	sprintf(buf, "%s.d", obj->savefile.name);
	relation_fp = fopen(buf, "wb");
	if (relation_fp == NULL) {
		logprintf(obj, "error: reldump can't open out file\n");
		exit(-1);
	}

	for (i = 0; i < filter->num_relations; i++) {
		fwrite(&r->rel_index, (size_t)1, sizeof(uint32), relation_fp);
		r = next_relation_ptr(r);
	}
	fclose(relation_fp);
}

/*--------------------------------------------------------------------*/
static void set_filtering_bounds(msieve_obj *obj, factor_base_t *fb, 
			uint32 filtmin_r, uint32 filtmin_a, 
			uint32 *entries_r_out, uint32 *entries_a_out,
			uint32 num_relations, uint32 force_small, 
			filter_t *filter) {

	uint32 entries_r, entries_a;

	if (force_small) {
		if (num_relations < 2000000)
			filtmin_r = filtmin_a = 30000;
		else if (num_relations < 10000000) {
			filtmin_r = MIN(filtmin_r / 2, 100000);
			filtmin_a = MIN(filtmin_a / 2, 100000);
		}
		else {
			filtmin_r = MIN(filtmin_r / 2, 720000);
			filtmin_a = MIN(filtmin_a / 2, 720000);
		}
	}

	logprintf(obj, "reading ideals above %u\n", filtmin_r);
	find_fb_size(fb, filtmin_r, filtmin_a, &entries_r, &entries_a);
	filter->filtmin_r = filtmin_r;
	filter->filtmin_a = filtmin_a;
	filter->target_excess = entries_r + entries_a;

	*entries_r_out = entries_r;
	*entries_a_out = entries_a;
}

/*--------------------------------------------------------------------*/
/* the multiple of the amount of excess needed for
   merging to proceed */

#define FINAL_EXCESS_FRACTION 1.16

static uint32 do_merge(msieve_obj *obj, filter_t *filter, merge_t *merge) {

	uint32 relations_needed;
	uint32 extra_needed = filter->target_excess;

	/* make the clique removal more conservative by leaving 
	   some of the excess; this makes the merge phase easier. 
	   Note that the singleton removal probably threw away 
	   many large ideals that occur too often to be worth 
	   tracking, which forces the target matrix size to 
	   increase, so that target_excess is larger now */

	filter->target_excess *= FINAL_EXCESS_FRACTION;

	if ((relations_needed = check_excess(filter)) > 0)
		return relations_needed;

	/* build the matrix */

	merge->num_extra_relations = NUM_EXTRA_RELATIONS;

	filter_make_relsets(obj, filter, merge, extra_needed);
	return 0;
}

/*--------------------------------------------------------------------*/
#define MAX_KEEP_WEIGHT 45

static uint32 do_partial_filtering(msieve_obj *obj, filter_t *filter,
				merge_t *merge, uint32 entries_r,
				uint32 entries_a) {

	uint32 relations_needed;
	uint32 max_weight = 20;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;

	if (filter->num_relations < 20000000)
		max_weight = 25;

	for (; max_weight < MAX_KEEP_WEIGHT; max_weight += 5) {

		filter->target_excess = entries_r + entries_a;
		filter->num_relations = num_relations;
		filter->num_ideals = num_ideals;

		filter_read_lp_file(obj, filter, max_weight);

		if ((relations_needed = do_merge(obj, filter, merge)) > 0)
			return relations_needed;

		/* retry the filtering with more aggressive
		   settings if the clique max weight was 
		   small and the largest number of relations 
		   in a relation set is close to the clique 
		   max weight. This is an indication that the 
		   clique processing ignored too many ideals 
		   that the merge phase should have seen */

		if (filter->max_ideal_weight < 18 &&
		    merge->max_relations > filter->max_ideal_weight - 3) {
			logprintf(obj, "matrix can improve, retrying\n");
			filter_free_relsets(merge);
			continue;
		}

		/* do not accept the collection of generated cycles 
		   unless the matrix they form is dense enough */

		if (merge->avg_cycle_weight > 63.0)
			break;

		logprintf(obj, "matrix not dense enough, retrying\n");
		filter_free_relsets(merge);
	}

	if (max_weight >= MAX_KEEP_WEIGHT) {
		printf("error: too many merge attempts\n");
		exit(-1);
	}
	return 0;
}

/*--------------------------------------------------------------------*/
uint32 nfs_filter_relations(msieve_obj *obj, mp_t *n) {

	filter_t filter;
	merge_t merge;
	uint32 filtmin_r, filtmin_a;
	uint32 entries_r, entries_a;
	uint32 num_relations;
	uint32 relations_needed;
	factor_base_t fb;
	time_t wall_time = time(NULL);
	uint64 savefile_size = get_file_size(obj->savefile.name);
	uint64 ram_size = (uint64)obj->mem_mb * 1048576;
	char lp_filename[256];

	logprintf(obj, "\n");
	logprintf(obj, "commencing relation filtering\n");

	if (ram_size == 0)
		ram_size = get_ram_size();

	memset(&filter, 0, sizeof(filter));
	memset(&merge, 0, sizeof(merge));
	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly, NULL)) {
		printf("filtering failed to read polynomials\n");
		exit(-1);
	}
	logprintf(obj, "estimated available RAM is %.1lf MB\n", 
				(double)ram_size / 1048576);

	/* delete duplicate relations */

	filtmin_r = filtmin_a = nfs_purge_duplicates(obj, &fb, 
					(uint32)obj->nfs_upper,
					&num_relations);

	/* set up the first disk-based pass; if the dataset is
	   "small", this will be the only such pass */

	set_filtering_bounds(obj, &fb, filtmin_r, filtmin_a,
				&entries_r, &entries_a, num_relations, 
				(uint32)(savefile_size < ram_size / 2), 
				&filter);

	/* separate out the large ideals and delete singletons
	   once they are all in memory. If the dataset is large,
	   first delete most of the singletons from the disk file */

	nfs_write_lp_file(obj, &fb, &filter, 
			(uint32)obj->nfs_upper, 0);
	if (filter.lp_file_size > ram_size / 2) {
		filter_purge_lp_singletons(obj, &filter, ram_size);
#if 0
		/* also delete most of the cliques from the disk
		   file if it is still large, and there is a
		   great deal of excess relations */

		if (filter.lp_file_size > ram_size / 2) {
			filter_purge_lp_cliques(obj, &filter);
		}
#endif
	}
	filter_read_lp_file(obj, &filter, 0);

	if (savefile_size < ram_size / 2) {

		/* dataset is "small"; build the matrix immediately. 
		   Depending on how much memory the machine has, really
		   big datasets may get to do this */

		if ((relations_needed = do_merge(obj, &filter, &merge)) > 0)
			return relations_needed;
	}
	else {  
		/* dataset is "large", perform multiple singleton passes.

		   The first pass used a large bound; the second 
		   filtering bound is much smaller. To allow reuse of 
		   previous results, the second bound is used 
		   during the rest of the filtering */

		dump_relation_numbers(obj, &filter);
		set_filtering_bounds(obj, &fb, filtmin_r, filtmin_a,
					&entries_r, &entries_a, 
					filter.num_relations, 1, &filter);

		free(filter.relation_array);
		filter.relation_array = NULL;

		nfs_write_lp_file(obj, &fb, &filter, 
				(uint32)obj->nfs_upper, 1);

		if (filter.lp_file_size < ram_size / 2) {

			/* dataset is small enough for filtering to 
			   complete in one pass */

			filter_read_lp_file(obj, &filter, 0);
			if ((relations_needed = do_merge(obj, &filter, 
							&merge)) > 0) {
				return relations_needed;
			}
		}
		else {
			/* dataset is so large that even the pruned version
			   cannot fit comfortably in memory. We have to put
			   up with reading only the ideals that occur in the
			   fewest relations, forming the matrix, and then
			   determining whether the matrix incorporates enough
			   of the dataset so that the matrix will work */

			if ((relations_needed = do_partial_filtering(obj,
						&filter, &merge, entries_r,
						entries_a)) > 0) {
				return relations_needed;
			}
		}
	}

	/* filtering succeeded; delete the LP file */

	sprintf(lp_filename, "%s.lp", obj->savefile.name);
	remove(lp_filename);

	/* optimize and then save the collection of relation-sets */

	filter_postproc_relsets(obj, &merge);
	filter_dump_relsets(obj, &merge);
	filter_free_relsets(&merge);
	wall_time = time(NULL) - wall_time;
	logprintf(obj, "RelProcTime: %u\n", (uint32)wall_time);
	return 0;
}
