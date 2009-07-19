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

#include "filter_priv.h"

#define MAX_WEIGHT 50

/*--------------------------------------------------------------------*/
void filter_purge_heavy_relations(msieve_obj *obj, 
				filter_t *filter) {

	/* combine deletion of the heaviest relations with 
	   in-memory singleton removal. We iterate until 
	   there are no more singletons and the difference between
	   the number of surviving relations and unique ideals
	   is close to filter->target_excess */

	uint32 i, j, k;
	uint32 *freqtable;
	uint32 weight_table[MAX_WEIGHT + 1] = {0};
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	relation_ideal_t *old_relation;
	uint32 orig_num_ideals;
	uint32 num_passes;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 new_num_relations;
	uint32 weight_target;
	uint32 delete_count;
	uint32 delete_target;

	logprintf(obj, "commencing heavy relation removal\n");

	num_relations = filter->num_relations;
	orig_num_ideals = num_ideals = filter->num_ideals;
	relation_array = filter->relation_array;
	freqtable = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* count the number of times each ideal occurs. Note
	   that since we know the exact number of ideals, we
	   don't need a hashtable to store the counts, just an
	   ordinary random-access array (i.e. a perfect hashtable) */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {

		uint32 curr_ideals = curr_relation->gf2_factors +
					curr_relation->ideal_count;

		weight_table[curr_ideals]++;

		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			freqtable[ideal]++;
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	logprintf(obj, "begin with %u relations and %u unique ideals\n", 
					num_relations, num_ideals);

	/* find the number of relations to delete and the
	   smallest weight of a deleted relation */

	delete_target = (filter->num_relations - filter->num_ideals) -
				filter->target_excess;
	if (delete_target > 1000)
		delete_target -= 1000;

	for (k = MAX_WEIGHT; k; k--) {
		if (weight_table[k])
			break;
	}
	weight_target = k;

	/* while relations keep getting deleted */

	delete_count =	num_passes = 0;
	new_num_relations = num_relations;
	do {
		num_relations = new_num_relations;
		new_num_relations = 0;
		curr_relation = relation_array;
		old_relation = relation_array;

		for (i = 0; i < num_relations; i++) {
			uint32 curr_num_ideals = curr_relation->ideal_count;
			uint32 delete = 0;
			uint32 ideal;
			relation_ideal_t *next_relation;

			next_relation = next_relation_ptr(curr_relation);

			k = curr_num_ideals + curr_relation->gf2_factors;

			if (k >= weight_target && 
			    delete_count < delete_target) {
				/* relation too heavy */
				delete = 1;
			}
			else {
				/* check the count of each ideal */

				for (j = 0; j < curr_num_ideals; j++) {
					ideal = curr_relation->ideal_list[j];
					if (freqtable[ideal] <= 1) {
						delete = 1;
						break;
					}
				}
			}

			if (delete) {

				/* decrement the count of each of the ideals.
				   If curr_relation contained the last
				   occurrence of a given ideal, then deleting
				   curr_relation will not change the target
				   excess */

				delete_count++;
				for (j = 0; j < curr_num_ideals; j++) {
					ideal = curr_relation->ideal_list[j];
					if (--freqtable[ideal] == 0)
						delete_target++;
				}

				/* if we've run out of relations with 
				   target_weight ideals, lower the weight 
				   that will be considered too heavy */

				if (--weight_table[k] == 0) {
					weight_target--;
					while (weight_target && weight_table[
							weight_target] == 0) {
						weight_target--;
					}
				}
			}
			else {
				/* relation survived this pass; append it to
				   the list of survivors */

				old_relation->rel_index = 
						curr_relation->rel_index;
				old_relation->gf2_factors = 
						curr_relation->gf2_factors;
				old_relation->ideal_count = curr_num_ideals;
				for (j = 0; j < curr_num_ideals; j++) {
					old_relation->ideal_list[j] =
						curr_relation->ideal_list[j];
				}
				new_num_relations++;
				old_relation = next_relation_ptr(old_relation);
			}

			curr_relation = next_relation;
		}

		num_passes++;
	} while (new_num_relations != num_relations);

	/* find the ideal that occurs in the most
	   relations, and renumber the ideals to ignore
	   any that have a count of zero */

	num_ideals = 0;
	for (i = j = 0; i < orig_num_ideals; i++) {
		if (freqtable[i]) {
			j = MAX(j, freqtable[i]);
			freqtable[i] = num_ideals++;
		}
	}

	logprintf(obj, "reduce to %u relations and %u ideals in %u passes\n", 
				num_relations, num_ideals, num_passes);
	logprintf(obj, "max relations containing the same ideal: %u\n", j);
	
	/* save the current state */

	filter->max_ideal_weight = j;
	filter->num_relations = num_relations;
	filter->num_ideals = num_ideals;
	filter->relation_array = relation_array = 
			(relation_ideal_t *)xrealloc(relation_array,
				(curr_relation - relation_array + 1) *
				sizeof(relation_ideal_t));

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			curr_relation->ideal_list[j] = freqtable[ideal];
		}
		curr_relation = next_relation_ptr(curr_relation);
	}
	free(freqtable);
}

/*--------------------------------------------------------------------*/
void filter_purge_singletons_core(msieve_obj *obj, 
				filter_t *filter) {

	/* main routine for performing in-memory singleton
	   removal. We iterate until there are no more singletons */

	uint32 i, j;
	uint32 *freqtable;
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	relation_ideal_t *old_relation;
	uint32 orig_num_ideals;
	uint32 num_passes;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 new_num_relations;

	logprintf(obj, "commencing in-memory singleton removal\n");

	num_relations = filter->num_relations;
	orig_num_ideals = num_ideals = filter->num_ideals;
	relation_array = filter->relation_array;
	freqtable = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* count the number of times each ideal occurs. Note
	   that since we know the exact number of ideals, we
	   don't need a hashtable to store the counts, just an
	   ordinary random-access array (i.e. a perfect hashtable) */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			freqtable[ideal]++;
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	logprintf(obj, "begin with %u relations and %u unique ideals\n", 
					num_relations, num_ideals);

	/* while singletons were found */

	num_passes = 0;
	new_num_relations = num_relations;
	do {
		num_relations = new_num_relations;
		new_num_relations = 0;
		curr_relation = relation_array;
		old_relation = relation_array;

		for (i = 0; i < num_relations; i++) {
			uint32 curr_num_ideals = curr_relation->ideal_count;
			uint32 ideal;
			relation_ideal_t *next_relation;

			next_relation = next_relation_ptr(curr_relation);

			/* check the count of each ideal */

			for (j = 0; j < curr_num_ideals; j++) {
				ideal = curr_relation->ideal_list[j];
				if (freqtable[ideal] <= 1)
					break;
			}

			if (j < curr_num_ideals) {

				/* relation is a singleton; decrement the
				   count of each of its ideals and skip it */

				for (j = 0; j < curr_num_ideals; j++) {
					ideal = curr_relation->ideal_list[j];
					freqtable[ideal]--;
				}
			}
			else {
				/* relation survived this pass; append it to
				   the list of survivors */

				old_relation->rel_index = 
						curr_relation->rel_index;
				old_relation->gf2_factors = 
						curr_relation->gf2_factors;
				old_relation->ideal_count = curr_num_ideals;
				for (j = 0; j < curr_num_ideals; j++) {
					old_relation->ideal_list[j] =
						curr_relation->ideal_list[j];
				}
				new_num_relations++;
				old_relation = next_relation_ptr(old_relation);
			}

			curr_relation = next_relation;
		}

		num_passes++;
	} while (new_num_relations != num_relations);

	/* find the ideal that occurs in the most
	   relations, and renumber the ideals to ignore
	   any that have a count of zero */

	num_ideals = 0;
	for (i = j = 0; i < orig_num_ideals; i++) {
		if (freqtable[i]) {
			j = MAX(j, freqtable[i]);
			freqtable[i] = num_ideals++;
		}
	}

	logprintf(obj, "reduce to %u relations and %u ideals in %u passes\n", 
				num_relations, num_ideals, num_passes);
	logprintf(obj, "max relations containing the same ideal: %u\n", j);
	
	/* save the current state */

	filter->max_ideal_weight = j;
	filter->num_relations = num_relations;
	filter->num_ideals = num_ideals;
	filter->relation_array = relation_array = 
			(relation_ideal_t *)xrealloc(relation_array,
				(curr_relation - relation_array + 1) *
				sizeof(relation_ideal_t));

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			curr_relation->ideal_list[j] = freqtable[ideal];
		}
		curr_relation = next_relation_ptr(curr_relation);
	}
	free(freqtable);
}
