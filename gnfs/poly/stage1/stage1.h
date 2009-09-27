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

#ifndef _STAGE1_H_
#define _STAGE1_H_

#include <poly_skew.h>

#ifdef __cplusplus
extern "C" {
#endif

#if MAX_POLY_DEGREE < 6
#error "supported poly degree must be at least 6"
#endif

#define MULTIPLIER 60	/* 2*2*3*5 */

/*-----------------------------------------------------------------------*/

/* search bounds */

typedef struct {
	mpz_t gmp_high_coeff_begin;
	mpz_t gmp_high_coeff_end;
	double norm_max; 
	double coeff_max;
	double p_size_max;
} bounds_t;

void stage1_bounds_init(bounds_t *bounds, poly_stage1_t *data);
void stage1_bounds_free(bounds_t *bounds);
void stage1_bounds_update(bounds_t *bounds, double N, 
			double high_coeff, uint32 degree);

/*-----------------------------------------------------------------------*/

typedef struct {

	uint32 degree;

	mpz_t high_coeff; 
	mpz_t N; 
	mpz_t trans_N;
	mpz_t trans_m0;
	mpz_t m0; 
	mpz_t p;
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_t tmp3;
	mpz_t tmp4;
	mpz_t tmp5;

	double coeff_max;
	double p_size_max;
	double sieve_size;

	stage1_callback_t callback;
	void *callback_data;
} poly_search_t;

void poly_search_init(poly_search_t *poly, poly_stage1_t *data);
void poly_search_free(poly_search_t *poly);

/*-----------------------------------------------------------------------*/

/* Rational leading coeffs of NFS polynomials are assumed 
   to be the product of two groups of factors p; each group 
   can be up to 64 bits in size and the product of (powers 
   of) up to MAX_P_FACTORS distinct primes */

#define MAX_P_FACTORS 4

#define MAX_ROOTS 8

/*-----------------------------------------------------------------------*/

typedef struct {
	uint32 p;
	uint32 r;
	uint8 log_p;
	uint8 num_roots;
	uint32 roots[MAX_POLY_DEGREE];
} sieve_prime_t;

typedef struct {
	sieve_prime_t *primes;
	uint32 num_primes;
	uint32 num_primes_alloc;
} sieve_prime_list_t;

typedef struct {

	uint32 degree;
	uint64 p_min, p_max;
	uint32 num_roots_min;
	uint32 num_roots_max;

	uint32 allow_prime_p;
	prime_sieve_t prime_sieve;
	uint64 next_prime_p;

	uint8 *sieve_block;
	uint32 curr_offset;
	sieve_prime_list_t good_primes;
	sieve_prime_list_t bad_primes;
	uint64 next_composite_p;

	mpz_t p, p2, m0, nmodp2, tmp1, tmp2;
	mpz_t accum[MAX_POLY_DEGREE + 1];
	mpz_t roots[MAX_ROOTS];
} sieve_fb_t;

void sieve_fb_init(sieve_fb_t *s, poly_search_t *poly,
			uint32 factor_min, uint32 factor_max);

void sieve_fb_free(sieve_fb_t *s);

void sieve_fb_reset(sieve_fb_t *s, uint64 p_min, uint64 p_max,
			uint32 num_roots_min, uint32 num_roots_max);

typedef void (*root_callback)(uint64 p, uint32 num_roots,
				mpz_t *roots, void *extra);

uint64 sieve_fb_next(sieve_fb_t *s, 
			poly_search_t *poly, 
			root_callback callback,
			void *extra);

/*-----------------------------------------------------------------------*/

/* main search routines */

void sieve_lattice(msieve_obj *obj, poly_search_t *poly, 
			uint32 small_fb_max, uint32 large_fb_min, 
			uint32 large_fb_max, uint32 deadline);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_H_ */
