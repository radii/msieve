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

#include "sqrt.h"

	/* This code computes the algebraic square root by
	   brute force. Given a collection of relations, each
	   of which is a polynomial p(x) = a*c-b*x for constant c
	   and various (a,b):

	   - Multiply all the relations together, modulo the
	     monic version of the algebraic polynomial f(x), 
	     of degree d. The result is a polynomial S(x) of degree
	     d-1 with arbitrary precision (very large) coefficients

	   - Find the polynomial that equals S(x) when squared
	     modulo f(x), using q-adic Newton iteration. This will
	     have coefficients about half the size of those in S(x)

	   - Convert the result to an integer modulo n by substituting
	     the root m of f(x) modulo n for x

	   None of the basic NFS papers consider this method to be
	   a viable option, but all of those papers were written
	   in the early to mid 1990s and computers were expensive and
	   limited then. On a modern machine, the memory consumed
	   by the direct method is not excessive (it's comparable to
	   the memory needed for the NFS linear algebra) and the 
	   runtime for optimized code using FFT-based multiplication
	   is quite reasonable. The big advantage of this method is
	   its simplicity; properly implementing the more standard
	   algorithms by Montgomery / Nguyen for the algebraic square
	   root would require approximately 100x as much code, plus
	   large amounts of very difficult algebraic number theory */
	     
/* for polynomials with arbitrary-precision coefficients.
   When the polynomial is monic, the highest-order coefficient
   (i.e. one) is considered implicit, and is *not* reflected
   in the degree */

typedef struct {
	uint32 degree;
	ap_t coeff[MAX_POLY_DEGREE];
} ap_poly_t;

/* bag of quantities needed for computing S(x) */

typedef struct {
	ap_poly_t *monic_poly;
	abpair_t *rlist;
	mp_t *c;
	fastmult_info_t *mult_info;
} relation_prod_t;

/*-------------------------------------------------------------------*/
static void ap_poly_init(ap_poly_t *poly) {
	
	uint32 i;

	for (i = poly->degree = 0; i < MAX_POLY_DEGREE; i++)
		ap_init(poly->coeff + i);
}

/*-------------------------------------------------------------------*/
static void ap_poly_clear(ap_poly_t *poly) {
	
	uint32 i;

	for (i = poly->degree = 0; i < MAX_POLY_DEGREE; i++)
		ap_clear(poly->coeff + i);
}

/*-------------------------------------------------------------------*/
static void ap_poly_mod_q(ap_poly_t *p, ap_t *q, ap_poly_t *res,
			ap_t *recip, fastmult_info_t *info) {

	uint32 i;

	for (i = 0; i <= p->degree; i++)
		ap_mod(p->coeff + i, q, recip, res->coeff + i, info);

	/* recalculate the degree */

	while (--i) {
		if (!ap_is_zero(&res->coeff[i]))
			break;
	}
	res->degree = i;
}

/*-------------------------------------------------------------------*/
static void ap_poly_bits(ap_poly_t *p, uint32 *total_bits, uint32 *max_bits) {

	uint32 i, bits1, bits2;

	for (i = bits1 = bits2 = 0; i <= p->degree; i++) {
		uint32 curr_bits = ap_bits(p->coeff + i);
		bits1 += curr_bits;
		bits2 = MAX(bits2, curr_bits);
	}

	*total_bits = bits1;
	*max_bits = bits2;
}

/*-------------------------------------------------------------------*/
static void ap_poly_monic_derivative(ap_poly_t *src, ap_poly_t *dest,
					fastmult_info_t *info) {
	
	uint32 i;
	ap_t c;

	/* compute the coefficients of the derivative of src,
	   assumed to be monic */

	ap_init(&c);
	for (i = 0; i < src->degree; i++) {
		ap_si2ap(i + 1, POSITIVE, &c);
		ap_mul(src->coeff + (i+1), &c, dest->coeff + i, info);
	}
	ap_si2ap(i + 1, POSITIVE, dest->coeff + i);

	dest->degree = src->degree;
	for (i++; i < MAX_POLY_DEGREE; i++)
		ap_clear(dest->coeff + i);
	ap_clear(&c);
}

/*-------------------------------------------------------------------*/
static void ap_poly_mul(ap_poly_t *p1, ap_poly_t *p2,
			ap_poly_t *mod, fastmult_info_t *info,
			uint32 free_p2) {

	/* multiply p1(x) by p2(x) modulo mod(x) (assumed monic)
	   If free_p2 is nonzero the coefficients of p2(x) are 
	   freed after being used */

	uint32 i, j;
	uint32 d = mod->degree;
	uint32 d1 = p1->degree;
	uint32 d2 = p2->degree;
	uint32 prod_degree;
	ap_t tmp[MAX_POLY_DEGREE + 1];
	ap_t prod;

	/* initialize */

	ap_init(&prod);
	for (i = 0; i < MAX_POLY_DEGREE + 1; i++)
		ap_init(tmp + i);

	/* multiply p1 by the leading coefficient of p2 */

	for (i = 0; i <= d1; i++) {
		ap_mul(p1->coeff + i, p2->coeff + d2, tmp + i, info);
	}
	prod_degree = d1;
	if (free_p2) {
		ap_clear(p2->coeff + d2);
	}

	/* for each of the other coefficients in p2 */

	for (i = d2 - 1; (int32)i >= 0; i--) {

		/* shift the accumulator up by one, zeroing
		   out the low-order coefficient. Recycle the
		   leading coefficent into the lowest-order
		   coefficient */

		ap_t high_coeff = tmp[prod_degree + 1];
		for (j = prod_degree; (int32)j >= 0; j--) {
			tmp[j + 1] = tmp[j];
		}
		tmp[0] = high_coeff;
		tmp[0].nwords = 0;

		/* add in the product of p1(x) and coefficient
		   i of p2 */

		for (j = d1; j; j--) {
			ap_mul(p1->coeff + j, p2->coeff + i, &prod, info);
			ap_add(tmp + j, &prod, tmp + j);
		}
		ap_mul(p1->coeff + j, p2->coeff + i, tmp + j, info);
		if (free_p2) {
			ap_clear(p2->coeff + i);
		}

		/* recalculate the degree of the result */

		prod_degree = d + 1;
		while (prod_degree && ap_is_zero(tmp + prod_degree))
			prod_degree--;

		/* if it exceeds the degree of mod(x), subtract
		   mod(x) * (leading accumulator coefficient) */

		if (prod_degree <= d)
			continue;

		for (j = d; (int32)j >= 0; j--) {
			ap_mul(mod->coeff + j, tmp + prod_degree, &prod, info);
			ap_sub(tmp + j, &prod, tmp + j);
		}
		prod_degree--;
	}

	/* move the result in the accumulator over to p1 */

	ap_clear(&prod);
	for (i = 0; i <= prod_degree; i++) {
		ap_clear(p1->coeff + i);
		p1->coeff[i] = tmp[i];
	}
	for (; i < MAX_POLY_DEGREE + 1; i++)
		ap_clear(tmp + i);

	/* recalculate the degree */

	i = prod_degree;
	while (i > 0 && ap_is_zero(p1->coeff + i)) {
		ap_clear(p1->coeff + i);
		i--;
	}
	p1->degree = i;

}

/*-------------------------------------------------------------------*/
static uint32 verify_product(ap_poly_t *ap_prod, abpair_t *abpairs, 
			uint32 num_relations, uint32 q, mp_t *c, 
			mp_poly_t *alg_poly) {

	/* a sanity check on the computed value of S(x): for
	   a small prime q for which alg_poly is irreducible,
	   verify that ap_prod mod q equals the product
	   mod q of the relations in abpairs[]. The latter can
	   be computed very quickly */

	uint32 i, j;
	uint32 c_mod_q = mp_mod_1(c, q);
	uint32 d = alg_poly->degree;
	uint32 prod[MAX_POLY_DEGREE];
	uint32 mod[MAX_POLY_DEGREE];
	uint32 accum[MAX_POLY_DEGREE + 1];
	ap_t ap_q;
	ap_poly_t tmp_poly;

	/* reduce ap_prod mod q */

	ap_init(&ap_q);
	ap_poly_init(&tmp_poly);
	ap_si2ap(q, POSITIVE, &ap_q);
	ap_poly_mod_q(ap_prod, &ap_q, &tmp_poly, NULL, NULL);

	/* compute the product mod q directly. First initialize
	   and reduce the coefficients of alg_poly mod q */

	for (i = 0; i < d; i++) {
		prod[i] = 0;
		mod[i] = mp_mod_1(&alg_poly->coeff[i].num, q);
		if (alg_poly->coeff[i].sign == NEGATIVE && mod[i] > 0) {
			mod[i] = q - mod[i];
		}
	}
	prod[0] = 1;

	/* multiply the product by each relation in
	   turn, modulo q */

	for (i = 0; i < num_relations; i++) {
		int64 a = abpairs[i].a;
		uint32 b = q - (abpairs[i].b % q);
		uint32 ac;

		a = a % (int64)q;
		if (a < 0)
			a += q;
		ac = mp_modmul_1((uint32)a, c_mod_q, q);

		for (j = accum[0] = 0; j < d; j++) {
			accum[j+1] = mp_modmul_1(prod[j], b, q);
			accum[j] = mp_modadd_1(accum[j],
					mp_modmul_1(ac, prod[j], q), q);
		}

		for (j = 0; j < d; j++) {
			prod[j] = mp_modsub_1(accum[j],
					mp_modmul_1(accum[d], mod[j], q), q);
		}
	}

	/* do the polynomial compare */

	for (i = 0; i < d; i++) {
		ap_t *coeff = tmp_poly.coeff + i;
		uint32 value = 0;
		if (!ap_is_zero(coeff)) {
			value = coeff->val[0];
			if (coeff->sign == NEGATIVE)
				value = q - value;
		}
		if (value != prod[i])
			break;
	}
	ap_poly_clear(&tmp_poly);
	ap_clear(&ap_q);
	if (i == d)
		return 1;
	return 0;
}

/*-------------------------------------------------------------------*/
static void relation_to_poly(abpair_t *abpair, mp_t *c,
				ap_poly_t *poly) {

	/* given a and b in abpair, along with c, compute
	   poly(x) = a*c - b*x */

	mp_t tmp;
	int64 abs_a;
	uint32 asign;
	mp_t mp_abs_a;

	abs_a = abpair->a;
	asign = POSITIVE;
	if (abs_a < 0) {
		abs_a = -abs_a;
		asign = NEGATIVE;
	}

	mp_abs_a.val[0] = (uint32)abs_a;
	mp_abs_a.val[1] = (uint32)(abs_a >> 32);
	mp_abs_a.nwords = 0;
	if (mp_abs_a.val[1])
		mp_abs_a.nwords = 2;
	else if (mp_abs_a.val[0])
		mp_abs_a.nwords = 1;

	mp_mul(c, &mp_abs_a, &tmp);
	ap_mp2ap(&tmp, asign, poly->coeff + 0);
	ap_si2ap(abpair->b, NEGATIVE, poly->coeff + 1);
	poly->degree = 1;
	if (ap_is_zero(poly->coeff + 1))
		poly->degree = 0;
}

/*-------------------------------------------------------------------*/
static void multiply_relations(relation_prod_t *prodinfo, 
			uint32 index1, uint32 index2,
			ap_poly_t *prod) {

	/* multiply together the relations from index1 
	   to index2, inclusive. We proceed recursively to
	   assure that polynomials with approximately equal-
	   size coefficients get multiplied, and also to
	   avoid wasting huge amounts of memory in the
	   beginning when all the polynomials are small
	   but the memory allocated for them is large */

	ap_poly_t prod1, prod2;

	if (index1 == index2) {
		/* base case of recursion */

		relation_to_poly(prodinfo->rlist + index1,
				 prodinfo->c, prod);
		return;
	}

	ap_poly_init(&prod1);
	ap_poly_init(&prod2);
	
	if (index1 == index2 - 1) {
		/* base case of recursion */

		relation_to_poly(prodinfo->rlist + index1,
				 prodinfo->c, &prod1);
		relation_to_poly(prodinfo->rlist + index2,
				 prodinfo->c, &prod2);
	}
	else {
		/* recursively compute the product of the first
		   half and the last half of the relations */

		uint32 mid = (index1 + index2) / 2;
		multiply_relations(prodinfo, index1, mid, &prod1);
		multiply_relations(prodinfo, mid + 1, index2, &prod2);
	}

	/* multiply them together and save the result */
	ap_poly_mul(&prod1, &prod2, prodinfo->monic_poly, 
			prodinfo->mult_info, 1);

	*prod = prod1;
}

/*-------------------------------------------------------------------*/
#define ISQRT_NUM_ATTEMPTS 10

static uint32 get_initial_inv_sqrt(msieve_obj *obj, mp_poly_t *mp_alg_poly,
				ap_poly_t *prod, ap_poly_t *isqrt_mod_q, 
				ap_t *q_out) {

	/* find the prime q_out and the initial value of the
	   reciprocal square root of prod(x) mod q_out to use 
	   for the Newton iteration */

	uint32 i, j;
	uint32 q, start_q;
	ap_poly_t prod_mod_q;
	mp_poly_t mp_prod_mod_q;
	mp_poly_t mp_isqrt_mod_q;

	/* find a prime q for which mp_alg_poly mod q is
	   irreducible. The starting value to try was passed in */

	ap_poly_init(&prod_mod_q);
	start_q = q_out->val[0];

	for (i = 0; i < ISQRT_NUM_ATTEMPTS; i++) {
		if (get_prime_for_sqrt(mp_alg_poly, start_q + 1, &q)) {
			logprintf(obj, "warning: no irreducible prime found, "
					"switching to small primes\n");
			if (start_q > 150)
				start_q = 50;
			get_prime_for_sqrt(mp_alg_poly, start_q + 1, &q);
		}

		ap_si2ap(q, POSITIVE, q_out);
		ap_poly_mod_q(prod, q_out, &prod_mod_q, NULL, NULL);

		/* convert prod_mod_q to an mp_poly_t */

		memset(&mp_prod_mod_q, 0, sizeof(mp_poly_t));
		mp_prod_mod_q.degree = prod_mod_q.degree;
		for (j = 0; j <= prod_mod_q.degree; j++) {
			signed_mp_t *coeff = mp_prod_mod_q.coeff + j;
			coeff->sign = prod_mod_q.coeff[j].sign;
			coeff->num.nwords = 1;
			coeff->num.val[0] = prod_mod_q.coeff[j].val[0];
		}

		/* find the reciprocal square root mod q, or try
		   another q if this fails */

		if (inv_sqrt_mod_q(&mp_isqrt_mod_q, &mp_prod_mod_q, 
				mp_alg_poly, q, &obj->seed1, &obj->seed2)) {
			break;
		}
		start_q = q;
	}

	ap_poly_clear(&prod_mod_q);
	if (i == ISQRT_NUM_ATTEMPTS) {
		logprintf(obj, "error: cannot recover square root mod q\n");
		return 0;
	}

	/* initialize isqrt_mod_q */

	for (i = 0; i <= mp_isqrt_mod_q.degree; i++) {
		signed_mp_t *coeff = mp_isqrt_mod_q.coeff + i;
		ap_mp2ap(&coeff->num, coeff->sign, isqrt_mod_q->coeff + i);
	}
	isqrt_mod_q->degree = mp_isqrt_mod_q.degree;

	logprintf(obj, "initial square root is modulo %u\n", q);
	return 1;
}

/*-------------------------------------------------------------------*/
#define DEFAULT_RBITS 20000000

static uint32 get_recip_bits(uint32 dbits, uint32 nbits,
				uint32 max_bits) {

	/* determine the number of bits in the reciprocal
	   used to compute num % den, assuming bit sizes
	   nbits and dbits, respectively. max_bits is the 
	   largest size of nbits to expect in the future. 
	   We can use this information to choose a reciprocal
	   size that avoids wasting time on tiny problems */

	double fft_size;
	uint32 fft_power, diff;
	
	/* if num is small or close to den in size, choose the 
	   reciprocal size to get the remainder in a single pass */

	if (dbits > nbits || nbits - dbits < DEFAULT_RBITS)
		return MIN(DEFAULT_RBITS, max_bits + 32);

	/* the more common case: num is larger than den, possibly
	   much larger. We don't want to make the reciprocal the 
	   size of num, since that can be huge. We can make it
	   sort of small, at the expense of lots of extra runtime
	   since remainders are calculated in several passes. A 
	   good compromise is to determine the smallest power-of-two
	   size FFT that will fit num, then choose the size of the
	   reciprocal so that the FFT is the same size.  */

	fft_size = (double)nbits / 16.0 + 4;
	fft_power = (uint32)(ceil(log(fft_size) / M_LN2));
	diff = (16 << fft_power) - nbits;
	if (diff > 128)
		diff -= 128;
	if (diff > nbits - dbits)
		diff = nbits - dbits;

	return MAX(diff, MIN(DEFAULT_RBITS, max_bits + 32));
}

/*-------------------------------------------------------------------*/
static uint32 get_final_sqrt(msieve_obj *obj, ap_poly_t *alg_poly,
			ap_poly_t *prod, ap_poly_t *isqrt_mod_q, 
			ap_t *q, fastmult_info_t *info) {

	/* the main q-adic Newton iteration. On input, isqrt_mod_q
	   contains the starting value of the reciprocal square
	   root R[0](x) of the polynomial prod(x). The iteration is

	   R[k](x) = R[k-1](x) * (3 - prod(x)*R[k-1](x)^2) / 2 mod (q^(2^k))

	   and at the end of iteration k, prod(x)*R[k-1](x)^2 mod (q^(2^k))
	   is 1. We keep iterating until q^(2^k) is larger than the
	   size of the coefficients of the square root (i.e. about half
	   the size of the coefficients of prod(x)). Then the square
	   root to use is R[k](x) * prod(x) mod (q^(2^k)), which is
	   written to isqrt_mod_q */

	uint32 i, j;
	uint32 qbits;
	uint32 rbits;
	uint32 prod_bits, prod_max_bits;
	uint32 num_iter;
	ap_poly_t tmp_poly;
	ap_t three;
	ap_t recip;
	uint32 converged = 0;

	/* initialize */

	ap_poly_init(&tmp_poly);
	ap_init(&three);
	ap_init(&recip);
	ap_si2ap(3, POSITIVE, &three);
	ap_poly_bits(prod, &prod_bits, &prod_max_bits);

	/* since prod(x) only matters mod q^(2^(final_k)), we can
	   cut the memory use in half by changing prod(x) to this.
	   Remember final_k as well */

	i = q->val[0];
	for (num_iter = 0; ap_bits(q) < prod_max_bits / 2 + 4000; num_iter++)
		ap_mul(q, q, q, info);

	rbits = get_recip_bits(ap_bits(q), prod_max_bits,
				prod_max_bits);
	ap_recip(q, &recip, rbits, info);
	ap_poly_mod_q(prod, q, prod, &recip, info);
	ap_si2ap(i, POSITIVE, q);
	ap_clear(&recip);

	/* trim the allocated memory for the new prod(x) */

	for (i = 0; i <= prod->degree; i++) {
		ap_t *coeff = prod->coeff + i;
		coeff->val = (uint32 *)xrealloc(coeff->val, 
					coeff->nwords * sizeof(uint32));
		coeff->num_alloc = coeff->nwords;
	}

	/* do the main iteration */

	for (i = 0; i < num_iter; i++) {

		/* square the previous modulus */

		ap_mul(q, q, q, info);
		qbits = ap_bits(q);

		/* if needed, compute a reciprocal to use in 
		   future mod operations */

		if (q->nwords > MAX_MP_WORDS) {
			rbits = get_recip_bits(qbits, 
						ap_bits(prod->coeff + 0),
						prod_max_bits);
			ap_recip(q, &recip, rbits, info);
		}

		/* compute prod(x) * (previous R)^2 */

		ap_poly_mod_q(prod, q, &tmp_poly, &recip, info);
		ap_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, info, 0);
		ap_poly_mod_q(&tmp_poly, q, &tmp_poly, &recip, info);
		ap_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, info, 0);
		ap_poly_mod_q(&tmp_poly, q, &tmp_poly, &recip, info);

		/* compute (3 - that) / 2 */

		ap_sub(&tmp_poly.coeff[0], &three, &tmp_poly.coeff[0]);
		for (j = 0; j <= tmp_poly.degree; j++) {
			ap_t *coeff = tmp_poly.coeff + j;

			if (coeff->sign == POSITIVE)
				coeff->sign = NEGATIVE;
			else
				coeff->sign = POSITIVE;

			/* thanks to Dario Alpern for pointing out how
			   easy division by 2 is when the modulus is odd */

			if (!ap_is_zero(coeff) && coeff->val[0] % 2)
				ap_add(coeff, q, coeff);
			ap_rshift(coeff, 1, coeff);
		}

		/* finally, compute the new R(x) by multiplying the
		   result above by the old R(x) */

		ap_poly_mul(&tmp_poly, isqrt_mod_q, alg_poly, info, 1);
		ap_poly_mod_q(&tmp_poly, q, isqrt_mod_q, &recip, info);
		ap_poly_clear(&tmp_poly);

		if (i == num_iter - 1) {
			/* this is the last iteration; attempt to compute
			   the square root. First multiply R(x) by prod(x),
			   deleting prod(x) since we won't need it beyond
			   this point */

			ap_poly_mul(isqrt_mod_q, prod, alg_poly, info, 1);
			ap_poly_mod_q(isqrt_mod_q, q, isqrt_mod_q, 
							&recip, info);

			/* this is a little tricky. Up until now we've
			   been working modulo big numbers, but the coef-
			   ficients of the square root are just integers,
			   and may be negative. Negative numbers mod q
			   have a numerical value near that of +q, but we
			   want the square root to have a negative coef-
			   ficient in that case. Hence, if the top
			   few words of any coefficent of the square root
			   match the top few words of q, we assume this
			   coefficient is negative and subtract q from it.

			   Theoretically we could be wrong, and the 
			   coefficient really is supposed to be a big 
			   positive number near q in size. However, if
			   q is thousands of bits larger than the size we
			   expect for the square root coefficients, this
			   is so unlikely that it's not worth worrying about */

			for (j = 0; j <= isqrt_mod_q->degree; j++) {
				ap_t *coeff = isqrt_mod_q->coeff + j;
				if (q->nwords == coeff->nwords &&
				    !memcmp(q->val + q->nwords - 4,
					    coeff->val + coeff->nwords - 4,
					    4 * sizeof(uint32))) {
					if (coeff->sign == NEGATIVE)
						ap_add(coeff, q, coeff);
					else
						ap_sub(coeff, q, coeff);
				}
			}

			/* another heuristic: we will assume the Newton
			   iteration has converged if, after applying the
			   correction above for negative square root
			   coefficients, the total number of bits in the 
			   coefficients of the resulting polynomial is
			   much smaller than we would expect from random
			   polynomials modulo q */

			ap_poly_bits(isqrt_mod_q, &prod_bits, &j);
			if (prod_bits < (isqrt_mod_q->degree + 1) * qbits - 100)
				converged = 1;
		}

		ap_clear(&recip);
	}

	/* clean up */

	if (!converged)
		logprintf(obj, "Newton iteration failed to converge\n");

	ap_poly_clear(&tmp_poly);
	ap_clear(&three);
	ap_clear(&recip);
	return converged;
}

/*-------------------------------------------------------------------*/
static void convert_to_integer(ap_poly_t *alg_sqrt, mp_t *n,
				mp_t *c, signed_mp_t *m1, 
				signed_mp_t *m0, mp_t *res) {

	/* given the completed square root, apply the homomorphism
	   to convert the polynomial to an integer. We do this
	   by evaluating alg_sqrt at c*m0/m1, with all calculations
	   performed mod n */

	uint32 i;
	ap_t ap_n;
	ap_t *ap_coeff;
	mp_t m1_pow;
	mp_t m1_tmp;
	mp_t m0_tmp;
	mp_t next_coeff;

	ap_init(&ap_n); 
	ap_mp2ap(n, POSITIVE, &ap_n);
	ap_poly_mod_q(alg_sqrt, &ap_n, alg_sqrt, NULL, NULL);
	ap_clear(&ap_n);

	mp_copy(&m1->num, &m1_tmp);
	if (m1->sign == NEGATIVE)
		mp_sub(n, &m1_tmp, &m1_tmp);
	mp_copy(&m1_tmp, &m1_pow);

	mp_modmul(&m0->num, c, n, &m0_tmp);
	if (m0->sign == POSITIVE)
		mp_sub(n, &m0_tmp, &m0_tmp);

	i = alg_sqrt->degree;
	ap_coeff = alg_sqrt->coeff + i;
	res->nwords = ap_coeff->nwords;
	memcpy(res->val, ap_coeff->val, res->nwords * sizeof(uint32));
	if (ap_coeff->sign == NEGATIVE)
		mp_sub(n, res, res);

	for (i--; (int32)i >= 0; i--) {
		mp_modmul(res, &m0_tmp, n, res);

		ap_coeff = alg_sqrt->coeff + i;
		mp_clear(&next_coeff);
		next_coeff.nwords = ap_coeff->nwords;
		memcpy(next_coeff.val, ap_coeff->val, 
				next_coeff.nwords * sizeof(uint32));
		mp_modmul(&next_coeff, &m1_pow, n, &next_coeff);
		if (ap_coeff->sign == NEGATIVE)
			mp_sub(n, &next_coeff, &next_coeff);

		mp_add(res, &next_coeff, res);
		if (i > 0)
			mp_modmul(&m1_pow, &m1_tmp, n, &m1_pow);
	}
	if (mp_cmp(res, n) > 0)
		mp_sub(res, n, res);
}

/*-------------------------------------------------------------------*/
void alg_square_root(msieve_obj *obj, mp_poly_t *mp_alg_poly, 
			mp_t *n, mp_t *c, signed_mp_t *m1, 
			signed_mp_t *m0, abpair_t *rlist, 
			uint32 num_relations, uint32 check_q,
			mp_t *sqrt_a) {
	
	/* external interface for computing the algebraic
	   square root */

	uint32 i;
	ap_poly_t alg_poly;
	ap_poly_t d_alg_poly;
	ap_poly_t prod;
	ap_poly_t alg_sqrt;
	fastmult_info_t info;
	relation_prod_t prodinfo;
	double log2_prodsize;
	ap_t *p;
	ap_t q;

	/* initialize */

	ap_init(&q);
	ap_poly_init(&alg_poly);
	ap_poly_init(&d_alg_poly);
	ap_poly_init(&prod);
	ap_poly_init(&alg_sqrt);
	fastmult_info_init(&info);

	/* convert the algebraic poly to arbitrary precision */

	for (i = 0; i < mp_alg_poly->degree; i++) {
		signed_mp_t *coeff = mp_alg_poly->coeff + i;
		ap_mp2ap(&coeff->num, coeff->sign, alg_poly.coeff + i);
	}
	alg_poly.degree = mp_alg_poly->degree - 1;

	/* multiply all the relations together */

	prodinfo.monic_poly = &alg_poly;
	prodinfo.rlist = rlist;
	prodinfo.c = c;
	prodinfo.mult_info = &info;

	logprintf(obj, "multiplying %u relations\n", num_relations);
	multiply_relations(&prodinfo, 0, num_relations - 1, &prod);
	logprintf(obj, "multiply complete, coefficients have about "
			"%3.2lf million bits\n",
			(double)ap_bits(prod.coeff + 0) / 1e6);

	/* perform a sanity check on the result */

	i = verify_product(&prod, rlist, num_relations, 
				check_q, c, mp_alg_poly);
	free(rlist);
	if (i == 0) {
		logprintf(obj, "error: relation product is incorrect\n");
		goto finished;
	}

	/* multiply by the square of the derivative of alg_poly;
	   this will guarantee that the square root of prod actually 
	   is an element of the number field defined by alg_poly.
	   If we didn't do this, we run the risk of the main Newton
	   iteration not converging */

	ap_poly_monic_derivative(&alg_poly, &d_alg_poly, &info);
	ap_poly_mul(&d_alg_poly, &d_alg_poly, &alg_poly, &info, 0);
	ap_poly_mul(&prod, &d_alg_poly, &alg_poly, &info, 1);

	/* pick the initial small prime to start the Newton iteration.
	   To save both time and memory, choose an initial prime 
	   such that squaring it a large number of times will produce
	   a value just a little larger than we need to calculate
	   the square root.
	
	   Note that contrary to what some authors write, pretty much
	   any starting prime is okay. The Newton iteration has a division
	   by 2, so that 2 must be invertible mod the prime (this is
	   guaranteed for odd primes). Also, the Newton iteration will
	   fail if both square roots have the same value mod the prime;
	   however, even a 16-bit prime makes this very unlikely */

	p = &prod.coeff[0];
	log2_prodsize = 32.0 * (p->nwords - 2) +
		        log(p->val[p->nwords - 1] * MP_RADIX +
			    p->val[p->nwords - 2]) / M_LN2 + 10000;
	while (log2_prodsize > 31.5)
		log2_prodsize *= 0.5;

	ap_si2ap((uint32)pow(2.0, log2_prodsize) + 1, POSITIVE, &q);

	/* get the initial inverse square root */

	if (!get_initial_inv_sqrt(obj, mp_alg_poly, 
				&prod, &alg_sqrt, &q)) {
		goto finished;
	}

	/* compute the actual square root */

	if (get_final_sqrt(obj, &alg_poly, &prod, &alg_sqrt, &q, &info))
		convert_to_integer(&alg_sqrt, n, c, m1, m0, sqrt_a);

finished:
	ap_poly_clear(&prod);
	ap_poly_clear(&alg_sqrt);
	ap_poly_clear(&alg_poly);
	ap_poly_clear(&d_alg_poly);
	ap_clear(&q);
	fastmult_info_free(&info);
}
