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

#include <common.h>
#include "gnfs.h"

/*--------------------------------------------------------------------*/
static uint32 divide_factor_out(mp_t *polyval, uint32 p, 
				uint32 *factors, uint32 *num_factors,
				uint32 compress) {

	/* read the rational factors. Note that the following
	   will work whether a given factor appears only once
	   or whether its full multiplicity is in the relation */

	uint32 i = *num_factors;
	uint32 multiplicity = 0;

	while (mp_mod_1(polyval, p) == 0) {
		mp_divrem_1(polyval, p, polyval);
		multiplicity++;
	}
	if (i + multiplicity >= TEMP_FACTOR_LIST_SIZE)
		return 1;

	if (compress) {
		if (multiplicity & 1)
			factors[i++] = p;
	}
	else if (multiplicity) {
		while (multiplicity--)
			factors[i++] = p;
	}
	*num_factors = i;
	return 0;
}

/*--------------------------------------------------------------------*/
int32 nfs_read_relation(char *buf, factor_base_t *fb, 
			relation_t *r, uint32 compress) {

	/* note that only the polynomials within the factor
	   base need to be initialized */

	uint32 i, p;
	int64 a, atmp;
	uint32 b;
	char *tmp, *next_field;
	signed_mp_t polyval;
	uint32 num_factors_r;
	uint32 num_factors_a;
	uint32 *factors = r->factors;

	/* read the relation coordinates */

	a = (int64)strtod(buf, &next_field);
	tmp = next_field;
	if (tmp[0] != ',' || !isdigit(tmp[1]))
		return -1;

	b = strtoul(tmp+1, &next_field, 10);
	tmp = next_field;

	num_factors_r = 0;
	num_factors_a = 0;
	r->a = a;
	r->b = b;

	/* for free relations, store the roots and not
	   the prime factors */

	if (b == 0) {
		uint32 i;
		uint32 roots[MAX_POLY_DEGREE];
		uint32 high_coeff, num_roots;

		p = (uint32)a;
		if (p == 0)
			return -2;

		num_roots = poly_get_zeros(roots, &fb->rfb.poly,
						p, &high_coeff, 0);
		if (num_roots != fb->rfb.poly.degree || high_coeff == 0)
			return -3;
		factors[num_factors_r++] = p;

		num_roots = poly_get_zeros(roots, &fb->afb.poly,
						p, &high_coeff, 0);
		if (num_roots != fb->afb.poly.degree || high_coeff == 0)
			return -4;
		for (i = 0; i < num_roots; i++)
			factors[num_factors_r + i] = roots[i];

		r->num_factors_r = num_factors_r;
		r->num_factors_a = num_roots;
		return 0;
	}

	if (tmp[0] != ':')
		return -5;
	
	atmp = a % (int64)b;
	if (atmp < 0)
		atmp += b;

	if (mp_gcd_1((uint32)atmp, b) != 1)
		return -6;

	/* handle a rational factor of -1 */

	eval_poly(&polyval, a, b, &fb->rfb.poly);
	if (mp_is_zero(&polyval.num))
		return -6;
	if (polyval.sign == NEGATIVE)
		factors[num_factors_r++] = 0;

	/* read the rational factors (possibly an empty list) */

	if (isxdigit(tmp[1])) {
		do {
			p = strtoul(tmp + 1, &next_field, 16);
			if (p > 1 && divide_factor_out(&polyval.num, p, 
						factors, 
						&num_factors_r, compress)) {
				return -8;
			}
			tmp = next_field;
		} while (tmp[0] == ',' && isxdigit(tmp[1]));
	}
	else {
		tmp++;
	}

	if (tmp[0] != ':')
		return -9;

	/* if there are rational factors still to be accounted
	   for, assume they are small and find them by brute force */

	for (i = p = 0; !mp_is_one(&polyval.num) && p < 1000; i++) {

		p += prime_delta[i];
		if (divide_factor_out(&polyval.num, p, factors,
					&num_factors_r, compress)) {
			return -10;
		}
	}

	if (!mp_is_one(&polyval.num))
		return -11;

	/* read the algebraic factors */

	eval_poly(&polyval, a, b, &fb->afb.poly);
	if (mp_is_zero(&polyval.num))
		return -12;

	if (isxdigit(tmp[1])) {
		do {
			p = strtoul(tmp + 1, &next_field, 16);
			if (p > 1 && divide_factor_out(&polyval.num, p, 
						factors + num_factors_r, 
						&num_factors_a, compress)) {
				return -13;
			}
			tmp = next_field;
		} while (tmp[0] == ',' && isxdigit(tmp[1]));
	}

	/* if there are algebraic factors still to be accounted
	   for, assume they are small and find them by brute force */

	for (i = p = 0; !mp_is_one(&polyval.num) && p < 1000; i++) {

		p += prime_delta[i];
		if (divide_factor_out(&polyval.num, p,
					factors + num_factors_r, 
					&num_factors_a, compress)) {
			return -14;
		}
	}

	if (!mp_is_one(&polyval.num))
		return -15;
	
	r->num_factors_r = num_factors_r;
	r->num_factors_a = num_factors_a;
	return 0;
}

/*--------------------------------------------------------------------*/
uint32 find_large_ideals(relation_t *rel, 
			relation_lp_t *out, 
			uint32 filtmin_r, uint32 filtmin_a) {
	uint32 i;
	uint32 num_ideals = 0;
	uint32 num_factors_r;
	int64 a = rel->a;
	uint32 b = rel->b;

	out->gf2_factors = 0;

	/* handle free relations */

	if (b == 0) {
		uint32 p = rel->factors[0];

		if (p > filtmin_r) {
			ideal_t *ideal = out->ideal_list + num_ideals;
			ideal->i.r = p;
			ideal->i.compressed_p = (p - 1) / 2;
			ideal->i.rat_or_alg = RATIONAL_IDEAL;
			num_ideals++;
		}
		else if (p > MAX_PACKED_PRIME) {
			out->gf2_factors++;
		}

		if (p > filtmin_a) {
			for (i = 0; i < rel->num_factors_a; i++) {
				ideal_t *ideal = out->ideal_list + 
							num_ideals + i;
				ideal->i.r = rel->factors[i + 1];
				ideal->i.compressed_p = (p - 1) / 2;
				ideal->i.rat_or_alg = ALGEBRAIC_IDEAL;
			}
			num_ideals += rel->num_factors_a;
		}
		else if (p > MAX_PACKED_PRIME) {
			out->gf2_factors += rel->num_factors_a;
		}

		out->ideal_count = num_ideals;
		return num_ideals;
	}

	/* find the large rational ideals */

	num_factors_r = rel->num_factors_r;

	for (i = 0; i < num_factors_r; i++) {
		uint32 p = rel->factors[i];

		/* if processing all the ideals, make up a
		   separate unique entry for rational factors of -1 */

		if (p == 0 && filtmin_r == 0) {
			ideal_t *ideal = out->ideal_list + num_ideals;
			ideal->i.compressed_p = 0x7fffffff;
			ideal->i.rat_or_alg = RATIONAL_IDEAL;
			ideal->i.r = (uint32)(-1);
			num_ideals++;
			continue;
		}

		if (p > filtmin_r) {

			/* make a single unique entry for p, instead
			   of finding the exact number r for which
			   rational_poly(r) mod p is zero */

			ideal_t *ideal = out->ideal_list + num_ideals;

			if (num_ideals >= TEMP_FACTOR_LIST_SIZE)
				return TEMP_FACTOR_LIST_SIZE + 1;

			ideal->i.compressed_p = (p - 1) / 2;
			ideal->i.rat_or_alg = RATIONAL_IDEAL;
			ideal->i.r = p;
			num_ideals++;
		}
		else if (p > MAX_PACKED_PRIME) {

			/* we only keep a count of the ideals that are
			   too small to list explicitly. NFS filtering
			   will work a little better if we completely
			   ignore the smallest ideals */

			out->gf2_factors++;
		}
	}

	/* repeat for the large algebraic ideals */

	for (i = 0; i < (uint32)rel->num_factors_a; i++) {
		uint32 p = rel->factors[num_factors_r + i];
		if (p > filtmin_a) {
			ideal_t *ideal = out->ideal_list + num_ideals;
			uint32 bmodp;

			if (num_ideals >= TEMP_FACTOR_LIST_SIZE)
				return TEMP_FACTOR_LIST_SIZE + 1;

			/* this time we have to find the exact r */

			bmodp = b % p;
			if (bmodp == 0) {
				ideal->i.r = p;
			}
			else {
				uint32 root;
				int64 mapped_a = a % (int64)p;
				if (mapped_a < 0)
					mapped_a += p;
				root = (uint32)mapped_a;
				root = mp_modmul_1(root, 
						mp_modinv_1(bmodp, p), p);
				ideal->i.r = root;
			}
			ideal->i.compressed_p = (p - 1) / 2;
			ideal->i.rat_or_alg = ALGEBRAIC_IDEAL;
			num_ideals++;
		}
		else if (p > MAX_PACKED_PRIME) {
			out->gf2_factors++;
		}
	}

	out->ideal_count = num_ideals;
	return num_ideals;
}

/*--------------------------------------------------------------------*/
static int compare_uint32(const void *x, const void *y) {
	uint32 *xx = (uint32 *)x;
	uint32 *yy = (uint32 *)y;
	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

static int bsearch_relation(const void *key, const void *rel) {
	relation_t *r = (relation_t *)rel;
	uint32 *k = (uint32 *)key;

	if ((*k) < r->rel_index)
		return -1;
	if ((*k) > r->rel_index)
		return 1;
	return 0;
}

static void nfs_get_cycle_relations(msieve_obj *obj, 
				factor_base_t *fb, uint32 num_cycles, 
				la_col_t *cycle_list, 
				uint32 *num_relations_out,
				relation_t **rlist_out,
				uint32 compress) {
	uint32 i, j;
	char buf[LINE_BUF_SIZE];
	relation_t *rlist;
	savefile_t *savefile = &obj->savefile;

	hashtable_t unique_relidx;
	uint32 num_unique_relidx;
	uint32 *relidx_list, *entry;

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	hashtable_init(&unique_relidx, (uint32)1, 0);

	/* fill the hashtable */

	for (i = 0; i < num_cycles; i++) {
		la_col_t *c = cycle_list + i;
		uint32 num_relations = c->cycle.num_relations;
		uint32 *list = c->cycle.list;

		for (j = 0; j < num_relations; j++) {
			hashtable_find(&unique_relidx, list + j, NULL, NULL);
		}
	}

	/* convert the internal list of hashtable entries into
	   a list of 32-bit relation numbers */

	hashtable_close(&unique_relidx);
	num_unique_relidx = hashtable_get_num(&unique_relidx);
	entry = (uint32 *)hashtable_get_first(&unique_relidx);
	relidx_list = unique_relidx.match_array;

	for (i = 0; i < num_unique_relidx; i++) {
		relidx_list[i] = *entry;
		entry = (uint32 *)hashtable_get_next(&unique_relidx, entry);
	}

	/* sort the list in order of increasing relation number */

	qsort(relidx_list, (size_t)num_unique_relidx, 
		sizeof(uint32), compare_uint32);

	logprintf(obj, "cycles contain %u unique relations\n", 
				num_unique_relidx);

	savefile_open(savefile, SAVEFILE_READ);

	/* read the list of relations */

	rlist = (relation_t *)xmalloc(num_unique_relidx * sizeof(relation_t));

	i = (uint32)(-1);
	j = 0;
	savefile_read_line(buf, sizeof(buf), savefile);
	while (!savefile_eof(savefile) && j < num_unique_relidx) {
		
		int32 status;

		if (buf[0] != '-' && !isdigit(buf[0])) {

			/* no relation at this line */

			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (++i < relidx_list[j]) {

			/* relation not needed */

			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		status = nfs_read_relation(buf, fb, &tmp_relation, compress);
		if (status) {
			/* at this point, if the relation couldn't be
			   read then the filtering stage should have
			   found that out and skipped it */

			logprintf(obj, "error: relation %u corrupt\n", i);
			exit(-1);
		}
		else {
			/* save the relation */

			uint32 num_r = tmp_relation.num_factors_r;
			uint32 num_a = tmp_relation.num_factors_a;
			relation_t *r = rlist + j++;

			*r = tmp_relation;
			r->rel_index = i;
			r->refcnt = 0;
			r->factors = (uint32 *)xmalloc((num_r + num_a) * 
							sizeof(uint32));
			memcpy(r->factors, tmp_relation.factors,
					(num_r + num_a) * sizeof(uint32));
		}

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	num_unique_relidx = *num_relations_out = j;
	logprintf(obj, "read %u relations\n", j);
	savefile_close(savefile);
	hashtable_free(&unique_relidx);

	/* walk through the list of cycles and convert
	   relation indices to relation pointers */

	for (i = 0; i < num_cycles; i++) {
		la_col_t *c = cycle_list + i;

		for (j = 0; j < c->cycle.num_relations; j++) {

			/* since relations were read in order of increasing
			   relation index (= savefile line number), use 
			   binary search to locate relation j for this
			   cycle, then save a pointer to it */

			relation_t *rptr = (relation_t *)bsearch(
						c->cycle.list + j,
						rlist,
						(size_t)num_unique_relidx,
						sizeof(relation_t),
						bsearch_relation);
			if (rptr == NULL) {
				/* this cycle is corrupt somehow */
				logprintf(obj, "error: cannot locate "
						"relation %u\n", 
						c->cycle.list[j]);
				exit(-1);
			}
			else {
				c->cycle.list[j] = rptr - rlist;
			}
		}
	}
	*rlist_out = rlist;
}

/*--------------------------------------------------------------------*/
void nfs_read_cycles(msieve_obj *obj, 
			factor_base_t *fb,
			uint32 *num_cycles_out, 
			la_col_t **cycle_list_out, 
			uint32 *num_relations_out,
			relation_t **rlist_out,
			uint32 compress,
			uint32 dependency) {

	uint32 i, j;
	uint32 num_cycles;
	uint32 num_relations;
	la_col_t *cycle_list;
	relation_t *rlist;

	/* read the raw list of relation numbers for each cycle */

	read_cycles(obj, &num_cycles, &cycle_list, dependency);

	if (num_cycles == 0) {
		free(cycle_list);
		*num_cycles_out = 0;
		*cycle_list_out = NULL;
		*num_relations_out = 0;
		*rlist_out = NULL;
		return;
	}

	/* give up if caller doesn't want the relations as well */

	if (fb == NULL || num_relations_out == NULL || rlist_out == NULL) {
		*num_cycles_out = num_cycles;
		*cycle_list_out = cycle_list;
		return;
	}

	/* now read the list of relations needed by the
	   list of cycles, and convert relation numbers
	   to relation pointers */

	nfs_get_cycle_relations(obj, fb, num_cycles, cycle_list, 
				&num_relations, &rlist, compress);

	/* count the number of times each relation is referenced */

	for (i = 0; i < num_cycles; i++) {
		la_col_t *c = cycle_list + i;
		for (j = 0; j < c->cycle.num_relations; j++) {
			relation_t *r = rlist + c->cycle.list[j];
			r->refcnt++;
		}
	}

	*num_cycles_out = num_cycles;
	*cycle_list_out = cycle_list;
	*num_relations_out = num_relations;
	*rlist_out = rlist;
}

/*--------------------------------------------------------------------*/
void nfs_free_relation_list(relation_t *rlist, uint32 num_relations) {

	uint32 i;

	for (i = 0; i < num_relations; i++)
		free(rlist[i].factors);
	free(rlist);
}
