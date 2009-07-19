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

/* Create <savefile_name>.s, a binary file containing the line
   numbers of relations that survived the singleton and clique
   removal phases.

   Removal of singletons happens in one of two modes, disk-based
   mode and in-memory mode. The disk-based mode uses a constant-size
   hashtable of bytes to store the counts of large ideals. Large
   ideals that hash to the same bin add their counts together, and
   singleton removal proceeds based on the hashtable until the number
   of singletons found decreases below a cutoff. This allows the first
   few singleton removal passes to run with constant (low) memory use, but
   some singletons may be missed.

   The in-memory mode runs either by itself or after the disk-based
   removal process completes. All the surviving large ideals get read
   in and assigned unique integers, and we use a perfect hash
   instead of an ordinary hash for the remaining singleton removal
   passes. In-memory singleton removal is rigorous, and using a perfect
   hashtable is much faster and more memory-efficient than the on-disk
   hashtable. The only drawback to in-memory mode is the memory use
   associated with the initial reading-in of large ideals */

/*--------------------------------------------------------------------*/
#define NUM_IDEAL_BINS 6

static uint8 *purge_singletons_pass1(msieve_obj *obj, factor_base_t *fb,
					filter_t *filter, 
					uint32 log2_hashtable_size,
					uint32 hash_mask, uint32 pass,
					uint32 *hash_bins_filled_out) {

	/* handle the initial phase of disk-based singleton
	   removal. This fills the hashtable initially */

	savefile_t *savefile = &obj->savefile;
	FILE *bad_relation_fp;
	uint32 i;
	char buf[LINE_BUF_SIZE];
	uint32 next_bad_relation;
	uint32 curr_relation;
	uint32 num_relations;
	uint8 *hashtable;
	uint32 hash_bins_filled;
	uint32 ideal_count_bins[NUM_IDEAL_BINS+2] = {0};

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass %u\n", pass);

	savefile_open(savefile, SAVEFILE_READ);

	if (pass == 1)
		sprintf(buf, "%s.d", savefile->name);
	else
		sprintf(buf, "%s.s", savefile->name);
	bad_relation_fp = fopen(buf, "rb");
	if (bad_relation_fp == NULL) {
		logprintf(obj, "error: singleton1 can't open rel file\n");
		exit(-1);
	}
	hashtable = (uint8 *)xcalloc(
			(size_t)1 << log2_hashtable_size, (size_t)1);

	/* for each relation declared good by the duplicate
	   removal phase */

	num_relations = 0;
	hash_bins_filled = 0;
	curr_relation = (uint32)(-1);
	next_bad_relation = (uint32)(-1);
	fread(&next_bad_relation, (size_t)1, 
			sizeof(uint32), bad_relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (++curr_relation == next_bad_relation) {
			fread(&next_bad_relation, (size_t)1, 
					sizeof(uint32), bad_relation_fp);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* tabulate the large ideals */

		num_ideals = find_large_ideals(&tmp_relation, &tmp_ideal, 
						filter->filtmin_r,
						filter->filtmin_a);

		/* increment the number of occurrences of each
		   large ideal in the hashtable. Each bin will
		   saturate at a count of 255 or more */

		for (i = 0; i < num_ideals; i++) {
			ideal_t *ideal = tmp_ideal.ideal_list + i;
			uint32 hashval = (HASH1(ideal->blob[0] ^ hash_mask) ^
					  HASH2(ideal->blob[1] ^ hash_mask)) >>
					(32 - log2_hashtable_size);

			if (hashtable[hashval] == 0)
				hash_bins_filled++;
			if (hashtable[hashval] < 0xff)
				hashtable[hashval]++;
		}

		/* update the histogram of ideal counts */

		if (num_ideals > NUM_IDEAL_BINS)
			ideal_count_bins[NUM_IDEAL_BINS+1]++;
		else
			ideal_count_bins[num_ideals]++;

		/* get the next line from the savefile */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* print some statistics */

	for (i = 0; i < NUM_IDEAL_BINS+1; i++) {
		logprintf(obj, "relations with %u large ideals: %u\n", 
					i, ideal_count_bins[i]);
	}
	logprintf(obj, "relations with %u+ large ideals: %u\n", 
				i, ideal_count_bins[i]);
	logprintf(obj, "%u relations and about %u large ideals\n", 
				num_relations, hash_bins_filled);

	savefile_close(savefile);
	fclose(bad_relation_fp);
	*hash_bins_filled_out = hash_bins_filled;
	return hashtable;
}

/*--------------------------------------------------------------------*/
static void purge_singletons_pass2(msieve_obj *obj, factor_base_t *fb,
				uint8 *hashtable, filter_t *filter,
				uint32 log2_hashtable_size, uint32 hash_mask,
				uint32 pass) {

	/* continue the disk-based singleton removal process */

	savefile_t *savefile = &obj->savefile;
	FILE *bad_relation_fp;
	FILE *out_fp;
	uint32 i;
	char buf[LINE_BUF_SIZE];
	char buf2[256];
	uint32 next_bad_relation;
	uint32 curr_relation;
	uint32 num_relations;
	uint32 num_singletons;
	uint32 num_ideals;

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass %u\n", pass);

	/* we start from either the file from the duplicate
	   removal or the file from the previous singleton
	   removal pass */

	savefile_open(savefile, SAVEFILE_READ);
	if (pass == 2)
		sprintf(buf, "%s.d", savefile->name);
	else
		sprintf(buf, "%s.s", savefile->name);
	bad_relation_fp = fopen(buf, "rb");
	if (bad_relation_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.s0", savefile->name);
	out_fp = fopen(buf, "wb");
	if (out_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open out file\n");
		exit(-1);
	}

	/* for each relation */

	num_relations = 0;
	num_singletons = 0;
	curr_relation = (uint32)(-1);
	next_bad_relation = (uint32)(-1);
	fread(&next_bad_relation, (size_t)1, 
			sizeof(uint32), bad_relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		uint32 is_singleton;
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (++curr_relation == next_bad_relation) {
			fwrite(&curr_relation, (size_t)1, 
					sizeof(uint32), out_fp);
			fread(&next_bad_relation, (size_t)1, 
					sizeof(uint32), bad_relation_fp);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* find the large ideals */

		num_ideals = find_large_ideals(&tmp_relation, 
					&tmp_ideal, 
					filter->filtmin_r,
					filter->filtmin_a);

		/* check the count of each large ideal */

		for (i = is_singleton = 0; i < num_ideals; i++) {
			ideal_t *ideal = tmp_ideal.ideal_list + i;
			uint32 hashval = (HASH1(ideal->blob[0] ^ hash_mask) ^
					  HASH2(ideal->blob[1] ^ hash_mask)) >>
					(32 - log2_hashtable_size);

			/* cache the hash values in case they 
			   are needed later */

			tmp_factors[i] = hashval;
			if (hashtable[hashval] <= 1)
				is_singleton = 1;
		}

		if (is_singleton) {

			/* decrement the frequency of all of the 
			   large ideals (unless a counter was 
			   saturated) */

			for (i = 0; i < num_ideals; i++) {
				uint32 hashval = tmp_factors[i];
				if (hashtable[hashval] < 0xff)
					hashtable[hashval]--;
			}
			num_singletons++;
			fwrite(&curr_relation, (size_t)1, 
				sizeof(uint32), out_fp);
		}

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* find the number of ideals remaining and
	   print some statistics */

	for (i = num_ideals = 0; 
			i < ((uint32)1 << log2_hashtable_size); i++) {
		if (hashtable[i] > 0)
			num_ideals++;
	}

	logprintf(obj, "found %u singletons\n", num_singletons);
	logprintf(obj, "current dataset: %u relations and "
			"about %u large ideals\n",
			num_relations - num_singletons, num_ideals);

	/* clean up and create the next singleton file */

	savefile_close(savefile);
	fclose(bad_relation_fp);
	fclose(out_fp);
	if (pass == 2) {
		sprintf(buf, "%s.d", savefile->name);
		remove(buf);
	}
	sprintf(buf, "%s.s0", savefile->name);
	sprintf(buf2, "%s.s", savefile->name);
	remove(buf2);
	if (rename(buf, buf2) != 0) {
		logprintf(obj, "error: singleton2 can't rename output file\n");
		exit(-1);
	}

	filter->num_relations = num_relations - num_singletons;
}

/*--------------------------------------------------------------------*/
static size_t read_lp_file(msieve_obj *obj, filter_t *filter, FILE *fp,
				size_t mem_use, 
				uint32 max_ideal_weight) {
	uint32 i, j, k;
	size_t header_size;
	relation_ideal_t tmp;
	relation_ideal_t *relation_array;
	uint32 curr_word;
	size_t num_relation_alloc;
	uint32 *counts;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;

	/* read in the relations from fp, as well as the large
	   ideals they contain. Do not save ideals that occur
	   more than max_ideal_weight times in the dataset */

	header_size = (sizeof(relation_ideal_t) - 
			sizeof(tmp.ideal_list)) / sizeof(uint32);
	counts = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* first build a frequency table for the large ideals */

	for (i = 0; i < num_relations; i++) {

		fread(&tmp, sizeof(uint32), header_size, fp);

		for (j = 0; j < tmp.ideal_count; j++) {
			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			counts[curr_ideal]++;
		}
	}

	/* renumber the ideals to ignore the ones that occur
	   too often */

	for (i = j = 0; i < num_ideals; i++) {
		if (counts[i] <= max_ideal_weight)
			counts[i] = j++;
		else
			counts[i] = (uint32)(-1);
	}
	filter->target_excess += i - j;
	filter->num_ideals = j;
	logprintf(obj, "keeping %u ideals with weight <= %u, "
			"new excess is %u\n",
			j, max_ideal_weight, 
			filter->target_excess);

	/* reread the relation list, saving the sparse ideals */

	rewind(fp);
	num_relation_alloc = 10000;
	curr_word = 0;
	relation_array = (relation_ideal_t *)xmalloc(
					num_relation_alloc *
					sizeof(relation_ideal_t));
	for (i = 0; i < num_relations; i++) {

		relation_ideal_t *r;

		/* make sure the relation array has room for the
		   new relation. Be careful increasing the array
		   size, since this is probably the largest array
		   in the NFS code */

		if (curr_word * sizeof(uint32) >=
				(num_relation_alloc-1) * 
				sizeof(relation_ideal_t)) {

			num_relation_alloc = 1.4 * num_relation_alloc;
			relation_array = (relation_ideal_t *)xrealloc(
					relation_array, 
					num_relation_alloc *
					sizeof(relation_ideal_t));
		}

		r = (relation_ideal_t *)(
			(uint32 *)relation_array + curr_word);
		fread(r, sizeof(uint32), header_size, fp);

		for (j = k = 0; j < r->ideal_count; j++) {

			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			curr_ideal = counts[curr_ideal];
			if (curr_ideal != (uint32)(-1))
				r->ideal_list[k++] = curr_ideal;
		}
		r->gf2_factors += j - k;
		r->ideal_count = k;
		curr_word += header_size + k;
	}

	/* finish up: trim the allocated relation array */

	filter->relation_array = (relation_ideal_t *)xrealloc(
						relation_array, 
						curr_word * 
						sizeof(uint32));
	free(counts);
	return MAX(mem_use, num_ideals * sizeof(uint32) +
			    num_relation_alloc * 
			    sizeof(relation_ideal_t));
}

/*--------------------------------------------------------------------*/
static void dump_relations(msieve_obj *obj, filter_t *filter) {

	uint32 i;
	char buf[256];
	FILE *relation_fp;
	relation_ideal_t *r = filter->relation_array;

	sprintf(buf, "%s.s", obj->savefile.name);
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
void nfs_purge_singletons(msieve_obj *obj, factor_base_t *fb,
			filter_t *filter, uint32 max_ideal_weight) {

	/* the last disk-based pass through the relation
	   file; its job is to form a packed array of 
	   relation_ideal_t structures */

	uint32 i;
	savefile_t *savefile = &obj->savefile;
	FILE *relation_fp;
	FILE *final_fp;
	char buf[LINE_BUF_SIZE];
	uint32 next_relation;
	uint32 curr_relation;
	uint32 num_relations;
	hashtable_t unique_ideals;
	size_t relation_array_words;
	size_t mem_use;
	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;
	uint32 have_bad_relation_list = (max_ideal_weight == 0);

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, final pass\n");

	savefile_open(savefile, SAVEFILE_READ);
	sprintf(buf, "%s.s", savefile->name);
	relation_fp = fopen(buf, "rb");
	if (relation_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.lp", savefile->name);
	final_fp = fopen(buf, "w+b");
	if (final_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open out file\n");
		exit(-1);
	}
	relation_array_words = 0;

	hashtable_init(&unique_ideals, (uint32)WORDS_IN(ideal_t), 0);

	/* for each relation that survived previous 
	   singleton filtering passes */

	curr_relation = (uint32)(-1);
	next_relation = (uint32)(-1);
	num_relations = 0;
	fread(&next_relation, (size_t)1, 
			sizeof(uint32), relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		int32 status;
		size_t curr_relation_words;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (have_bad_relation_list) {
			if (++curr_relation == next_relation) {
				fread(&next_relation, (size_t)1, 
						sizeof(uint32), relation_fp);
				savefile_read_line(buf, sizeof(buf), savefile);
				continue;
			}
		}
		else {
			if (++curr_relation < next_relation) {
				savefile_read_line(buf, sizeof(buf), savefile);
				continue;
			}
			fread(&next_relation, (size_t)1, 
					sizeof(uint32), relation_fp);
		}

		/* read it in */

		status = nfs_read_relation(buf, fb, &tmp_relation, 1);

		if (status == 0) {
			relation_lp_t tmp_ideal;
			relation_ideal_t packed_ideal;
			num_relations++;

			/* get the large ideals */

			find_large_ideals(&tmp_relation, &tmp_ideal, 
						filter->filtmin_r,
						filter->filtmin_a);

			packed_ideal.rel_index = curr_relation;
			packed_ideal.gf2_factors = tmp_ideal.gf2_factors;
			packed_ideal.ideal_count = tmp_ideal.ideal_count;

			/* map each ideal to a unique integer */

			for (i = 0; i < tmp_ideal.ideal_count; i++) {
				ideal_t *ideal = tmp_ideal.ideal_list + i;

				hashtable_find(&unique_ideals, ideal,
						packed_ideal.ideal_list + i,
						NULL);
			}

			/* dump the relation to disk */

			curr_relation_words = (sizeof(packed_ideal) -
					(TEMP_FACTOR_LIST_SIZE - 
					 	packed_ideal.ideal_count) * 
					sizeof(packed_ideal.ideal_list[0])) /
					sizeof(uint32);
			fwrite(&packed_ideal, curr_relation_words,
				sizeof(uint32), final_fp);
			relation_array_words += curr_relation_words;
		}

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	savefile_close(savefile);
	fclose(relation_fp);

	/* finish up; the disk-based pass simply reads all
	   the large ideals from disk, but subsequent passes
	   filter out large ideals that occur too often */

	filter->num_relations = num_relations;
	filter->num_ideals = hashtable_get_num(&unique_ideals);
	mem_use = hashtable_sizeof(&unique_ideals);
	hashtable_free(&unique_ideals);

	rewind(final_fp);
	if (max_ideal_weight == 0) {
		filter->relation_array = (relation_ideal_t *)xmalloc(
						relation_array_words * 
						sizeof(uint32));
		fread(filter->relation_array, relation_array_words,
				sizeof(uint32), final_fp);
		mem_use = MAX(mem_use, 
				relation_array_words * sizeof(uint32));
	}
	else {
		mem_use = read_lp_file(obj, filter, final_fp, 
					mem_use, max_ideal_weight);
	}

	logprintf(obj, "memory use: %.1f MB\n",
			(double)mem_use / 1048576);
	fclose(final_fp);
	sprintf(buf, "%s.lp", savefile->name);
	remove(buf);

	/* remove the rest of the singletons */

	filter_purge_singletons_core(obj, filter);
}

/*--------------------------------------------------------------------*/
void nfs_purge_singletons_initial(msieve_obj *obj, 
				factor_base_t *fb, filter_t *filter) {

	uint32 hash_mask = 0x0f0f0f0f;
	uint32 log2_hashtable_size = 27;
	uint32 hash_bins_filled = 0;
	uint32 curr_pass = 1;
	uint32 num_relations;
	uint8 *hashtable = purge_singletons_pass1(obj, fb, filter, 
						log2_hashtable_size,
						hash_mask, curr_pass++, 
						&hash_bins_filled);

	/* perform disk-based singleton removal until the
	   number of relations is small or the number of 
	   singletons found becomes small */

	filter->num_relations = 0xffffffff;
	do {
		num_relations = filter->num_relations;
		purge_singletons_pass2(obj, fb, hashtable, 
					filter, log2_hashtable_size,
					hash_mask, curr_pass++);

		/* if the hashtable was saturated going 
		   into the above singleton removal pass, 
		   rebuild it with more entries and change
		   the hash function */

		if (hash_bins_filled > 
			0.3 * (1 << log2_hashtable_size)) {
			free(hashtable);
			log2_hashtable_size = MIN(30, log2_hashtable_size + 1);
			hash_mask = (hash_mask << 1) | (hash_mask >> 31);
			hashtable = purge_singletons_pass1(
						obj, fb, filter, 
						log2_hashtable_size,
						hash_mask, curr_pass++, 
						&hash_bins_filled);
		}
	} while (filter->num_relations > 200000 &&
		num_relations - filter->num_relations > 500000);

	free(hashtable);

	nfs_purge_singletons(obj, fb, filter, 0);

	/* save the relation list once it is free of singletons.
	   This will convert the list of relations-to-skip into a list
	   of relations-to-keep */

	dump_relations(obj, filter);
}

