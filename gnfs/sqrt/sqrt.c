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

#include "sqrt.h"

/* we will need to find primes q for which f(x) mod q
   is irreducible. This is the maximum number of q to try */

#define NUM_PRIME_RETRIES 100

/*--------------------------------------------------------------------*/
uint32 get_prime_for_sqrt(mp_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out) {

	uint32 i;
	uint32 status = 0;
	uint32 q = 0;
	prime_sieve_t prime_sieve;

	init_prime_sieve(&prime_sieve, min_value, (uint32)(-1));
	for (i = 0; i < NUM_PRIME_RETRIES; i++) {
		uint32 tmp_q = get_next_prime(&prime_sieve);
		if (is_irreducible(alg_poly, tmp_q)) {
			q = tmp_q;
			break;
		}
		else if (q == 0) {
			/* in rare cases, there is no q for which alg_poly
			   mod q is irreducible. Life becomes much more
			   difficult in this case, since square roots mod
			   q will not be unique. The only alternative when
			   this happens is to pick a q such that alg_poly 
			   mod q has no linear roots (so that all of the
			   relations mod q are relatively prime to alg_poly), 
			   then keep trying dependencies until by luck we 
			   start with the correct initial square root for 
			   the Newton iteration to succeed.

			   Buhler et. al. show that for polynomial degree d,
			   on average one in 2^(d/2) dependencies will lead to
			   a congruence of squares (and about half of those
			   will lead to a factor). Technically we also need to
			   check that alg_poly mod q is squarefree, but that
			   would require a full polynomial factoring routine;
			   I'm gambling that being squarefree is not rare. */

			uint32 roots[MAX_POLY_DEGREE];
			uint32 high_coeff;
			uint32 num_roots = poly_get_zeros(roots, alg_poly,
						tmp_q, &high_coeff, 1);
			if (high_coeff != 0 && num_roots == 0)
				q = tmp_q;
		}
	}
	free_prime_sieve(&prime_sieve);

	if (i == NUM_PRIME_RETRIES)
		status = 1;

	*q_out = q;
	return status;
}

/*--------------------------------------------------------------------*/
static void eval_poly_derivative(mp_poly_t *poly, 
				signed_mp_t *m1, signed_mp_t *m0, 
				mp_t *n, mp_t *res) {
	uint32 i;
	mp_t next_coeff;
	mp_t m1_pow, m1_tmp, m0_tmp;

	mp_copy(&m1->num, &m1_tmp);
	if (m1->sign == NEGATIVE)
		mp_sub(n, &m1_tmp, &m1_tmp);
	mp_copy(&m1_tmp, &m1_pow);

	mp_copy(&m0->num, &m0_tmp);
	if (m0->sign == POSITIVE)
		mp_sub(n, &m0_tmp, &m0_tmp);

	i = poly->degree;
	mp_mul_1(&poly->coeff[i].num, i, res);

	while (--i) {
		signed_mp_t *coeff = poly->coeff + i;

		mp_modmul(res, &m0_tmp, n, res);

		mp_mul_1(&coeff->num, i, &next_coeff);
		if (coeff->sign == NEGATIVE)
			mp_sub(n, &next_coeff, &next_coeff);
		mp_modmul(&next_coeff, &m1_pow, n, &next_coeff);
		mp_add(res, &next_coeff, res);

		if (i > 1)
			mp_modmul(&m1_pow, &m1_tmp, n, &m1_pow);
	}

	if (mp_cmp(res, n) >= 0)
		mp_sub(res, n, res);
}

/*--------------------------------------------------------------------*/
typedef struct {
	uint32 p;
	uint32 count;
} rat_prime_t;

static uint32 rat_square_root(relation_t *rlist, uint32 num_relations,
				mp_t *n, mp_t *sqrt_r) {
	uint32 i, j, num_primes;
	hashtable_t h;
	uint32 already_seen;
	uint32 array_size;
	mp_t base, exponent, tmp;
	uint32 status = 0;
	rat_prime_t *curr;

	/* count up the number of times each prime factor in
	   rlist occurs */

	hashtable_init(&h, (uint32)WORDS_IN(rat_prime_t), 1);

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;

		for (j = array_size = 0; j < r->num_factors_r; j++) {
			uint32 p = decompress_p(r->factors, &array_size);
			curr = (rat_prime_t *)hashtable_find(&h, &p, NULL,
							    &already_seen);
			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}
	}

	/* verify all such counts are even, and form the 
	   rational square root */

	mp_clear(&base); base.nwords = 1;
	mp_clear(&exponent); exponent.nwords = 1;
	mp_clear(sqrt_r); sqrt_r->nwords = sqrt_r->val[0] = 1;
	num_primes = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_primes; i++) {
		if (curr->count % 2) {
			status = 1;
			break;
		}
		if (curr->p && curr->count) {
			base.val[0] = curr->p;
			exponent.val[0] = curr->count / 2;
			mp_expo(&base, &exponent, n, &tmp);
			mp_modmul(sqrt_r, &tmp, n, sqrt_r);
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	return status;
}

/*--------------------------------------------------------------------*/
typedef struct {
	ideal_t ideal;
	uint32 count;
} alg_prime_t;

static uint32 verify_alg_ideal_powers(relation_t *rlist, 
					uint32 num_relations,
					uint32 *num_free_relations) {

	uint32 i, j, num_ideals;
	hashtable_t h;
	uint32 already_seen;
	alg_prime_t *curr;
	uint32 status = 0;

	/* count the multiplicity of each algebraic ideal (not
	   just the prime to which the ideal corresponds) in rlist */

	*num_free_relations = 0;

	hashtable_init(&h, (uint32)WORDS_IN(alg_prime_t),
			(uint32)WORDS_IN(ideal_t));

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		relation_lp_t rlp;

		find_large_ideals(r, &rlp, 0, 0);

		for (j = 0; j < rlp.ideal_count; j++) {
			ideal_t *curr_ideal = rlp.ideal_list + j;

			if (curr_ideal->i.rat_or_alg == RATIONAL_IDEAL)
				continue;

			curr = (alg_prime_t *)hashtable_find(&h, curr_ideal, 
						NULL, &already_seen);

			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}

		if (r->b == 0)
			(*num_free_relations)++;
	}

	/* verify each ideal occurs an even number of times */

	num_ideals = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_ideals; i++) {
		if (curr->count % 2) {
			status = 1;
			break;
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	return status;
}

/*--------------------------------------------------------------------*/
uint32 nfs_find_factors(msieve_obj *obj, mp_t *n, 
			factor_list_t *factor_list) {

	/* external interface for the NFS square root */

	uint32 i, j;
	uint32 check_q;
	factor_base_t fb;
	mp_poly_t monic_alg_poly;
	mp_t sqrt_r, sqrt_a;
	mp_t c, tmp1, tmp2;
	signed_mp_t *m0;
	signed_mp_t *m1;
	uint32 dep_lower = 1;
	uint32 dep_upper = 64;
	uint32 factor_found = 0;
	time_t cpu_time;

	logprintf(obj, "\n");
	logprintf(obj, "commencing square root phase\n");

	/* read in the NFS polynomials */

	cpu_time = time(NULL);
	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly, NULL)) {
		logprintf(obj, "polynomials not found\n");
		return 0;
	}

	/* find the values needed to convert the algebraic 
	   square root back to an integer */

	if (fb.rfb.poly.degree != 1) {
		logprintf(obj, "cannot handle non-linear polynomials\n");
		return 0;
	}
	m0 = &fb.rfb.poly.coeff[0];
	m1 = &fb.rfb.poly.coeff[1];

	/* construct a monic version of the algebraic poly,
	   saving off the leading coefficient separately */

	j = fb.afb.poly.degree;
	if (fb.afb.poly.coeff[j].sign == NEGATIVE) {
		logprintf(obj, "cannot handle negative leading "
				"algebraic polynomial coefficient\n");
		return 0;
	}
	memset(&monic_alg_poly, 0, sizeof(mp_poly_t));
	mp_copy(&fb.afb.poly.coeff[j].num, &c);
	mp_copy(&c, &tmp1);
	monic_alg_poly.coeff[j-1] = fb.afb.poly.coeff[j-1];
	for (i = j - 2; (int32)i >= 0; i--) {
		mp_mul(&fb.afb.poly.coeff[i].num, &tmp1,
				&monic_alg_poly.coeff[i].num);
		monic_alg_poly.coeff[i].sign = fb.afb.poly.coeff[i].sign;

		if (i > 0) {
			mp_mul(&c, &tmp1, &tmp2);
			mp_copy(&tmp2, &tmp1);
		}
	}
	monic_alg_poly.coeff[j].num.nwords = 1;
	monic_alg_poly.coeff[j].num.val[0] = 1;
	monic_alg_poly.degree = fb.afb.poly.degree;
	get_prime_for_sqrt(&monic_alg_poly, (uint32)0x80000000, &check_q);

	/* determine the list of dependencies to compute */

	if (obj->nfs_lower && obj->nfs_upper) {
		dep_lower = MIN(obj->nfs_lower, 64);
		dep_upper = MIN(obj->nfs_upper, 64);
		dep_upper = MAX(dep_lower, dep_upper);
	}

	/* for each dependency */

	for (i = dep_lower; i <= dep_upper; i++) {

		uint32 num_relations;
		uint32 num_free_relations;
		relation_t *rlist;
		abpair_t *abpairs;
		mp_t exponent;

		logprintf(obj, "reading relations for dependency %u\n", i);

		/* read in only the relations for dependency i */

		nfs_read_cycles(obj, &fb, NULL, NULL,
				&num_relations, &rlist, 0, i);

		if (num_relations == 0)
			continue;

		/* do some sanity checking, performing increasing
		   amounts of work as the dependency proves itself
		   to be valid */

		if (num_relations % 2) {
			/* the number of relations in the dependency must
			   be even, because each relation represents a
			   degree-1 polynomial, and the product of these
			   relations will not have a square root unless the
			   degree of the product is even */
			logprintf(obj, "number of relations is not even\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (verify_alg_ideal_powers(rlist, 
				num_relations, &num_free_relations) != 0) {
			logprintf(obj, "algebraic side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (num_free_relations % 2) {
			logprintf(obj, "number of free relations (%u) is "
					"not even\n", num_free_relations);
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (rat_square_root(rlist, num_relations, n, &sqrt_r) != 0) {
			logprintf(obj, "rational side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}

		/* flatten the list of relations; each occurrence of
		   a relation gets its own abpair_t */

		abpairs = (abpair_t *)xmalloc(num_relations *
						sizeof(abpair_t));
		for (j = 0; j < num_relations; j++) {
			abpairs[j].a = rlist[j].a;
			abpairs[j].b = rlist[j].b;
		}
		nfs_free_relation_list(rlist, num_relations);

		/* perform the major work: the algebraic square root.
		   Note that to conserve memory, abpairs is freed in
		   the following call */

		mp_clear(&sqrt_a);
		alg_square_root(obj, &monic_alg_poly, n, &c, m1, m0, abpairs, 
					num_relations, check_q, &sqrt_a);
		if (mp_is_zero(&sqrt_a)) {
			logprintf(obj, "algebraic square root failed\n");
			continue;
		}

		/* an algebraic square root is available; move on
		   to the final congruence of squares. The arithmetic
		   is as given in Buhler et. al. with one exception:
		   when the rational poly is nonmonic there is a 
		   correction to the final square root value but the 
		   free relations *do not* figure into it. This latter
		   point is completely ignored in the literature! */

		eval_poly_derivative(&fb.afb.poly, m1, m0, n, &tmp1);
		mp_modmul(&tmp1, &sqrt_r, n, &sqrt_r);

		mp_clear(&exponent);
		exponent.nwords = 1;
		if (!mp_is_one(&c)) {
			exponent.val[0] = num_relations/ 2 + 
						fb.afb.poly.degree - 2;
			mp_expo(&c, &exponent, n, &tmp2);
			mp_modmul(&tmp2, &sqrt_r, n, &sqrt_r);
		}

		if (!mp_is_one(&m1->num)) {
			exponent.val[0] = (num_relations -
				       		num_free_relations) / 2;
			mp_copy(&m1->num, &tmp1);
			if (m1->sign == NEGATIVE)
				mp_sub(n, &tmp1, &tmp1);
			mp_expo(&tmp1, &exponent, n, &tmp2);
			mp_modmul(&tmp2, &sqrt_a, n, &sqrt_a);
		}

		/* a final sanity check: square the rational and algebraic 
		   square roots, expecting the same value modulo n */

		mp_modmul(&sqrt_r, &sqrt_r, n, &tmp1);
		mp_modmul(&sqrt_a, &sqrt_a, n, &tmp2);
		if (mp_cmp(&tmp1, &tmp2) != 0) {
			logprintf(obj, "dependency does not form a "
					"congruence of squares!\n");
			continue;
		}

		/* look for a nontrivial factor of n */

		mp_add(&sqrt_r, &sqrt_a, &tmp1);
		mp_gcd(&tmp1, n, &tmp1);
		if (!mp_is_one(&tmp1) && mp_cmp(n, &tmp1) != 0) {
			/* factor found; add it to the list of factors. 
			   Stop trying dependencies if the remaining
			   composite is small enough that another method
			   will factor it faster.

			   Actually, we should be stopping when the remaining
			   composite is much larger (70-80 digits), but 
			   avoid doing this because the MPQS code will run
			   and wipe out all the NFS relations we've collected */

			factor_found = 1;
			if (factor_list_add(obj, factor_list, &tmp1) < 
						SMALL_COMPOSITE_CUTOFF_BITS)
				break;
		}
	}

	cpu_time = time(NULL) - cpu_time;
	logprintf(obj, "sqrtTime: %u\n", (uint32)cpu_time);
	return factor_found;
}
