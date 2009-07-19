#ifndef _STAGE1_H_
#define _STAGE1_H_

#include <poly_skew.h>

#ifdef __cplusplus
extern "C" {
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

/* data for one polynomial */

typedef struct {
	mpz_t high_coeff; 
	mpz_t trans_N;
	mpz_t trans_m0;
	double coeff_max;
	double p_size_max;

	double sieve_size;
} curr_poly_t;

/* data for a batch of polynomials */

#define POLY_BATCH_SIZE 50

typedef struct {
	uint32 num_poly;
	curr_poly_t batch[POLY_BATCH_SIZE];

	mpz_t N; 
	mpz_t m0; 
	mpz_t p;
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_t tmp3;
	mpz_t tmp4;
	mpz_t tmp5;

	uint32 degree;
	stage1_callback_t callback;
	void *callback_data;
} poly_batch_t;

void poly_batch_init(poly_batch_t *poly, poly_stage1_t *data);
void poly_batch_free(poly_batch_t *poly);

/*-----------------------------------------------------------------------*/
typedef struct {
	uint16 which_poly;
	uint32 start_offset;
} index_t;

typedef struct {
	uint32 p;
	uint16 num_poly;
	uint16 num_roots;
	uint32 index_start_offset;
} p_entry_t;

typedef struct {
	uint32 num_entries;
	uint32 num_small_entries;
	uint32 num_entries_alloc;
	p_entry_t *entries;

	uint32 num_index;
	uint32 num_index_alloc;
	index_t *index;

	uint32 num_roots;
	uint32 num_roots_alloc;
	uint64 *roots;
} p_batch_t;

void p_batch_init(p_batch_t *pbatch);
void p_batch_free(p_batch_t *pbatch);
void p_batch_reset(p_batch_t *pbatch);
void p_batch_add(p_batch_t *pbatch, uint64 p, uint32 poly_idx,
			uint32 num_roots, mpz_t *roots);

/*-----------------------------------------------------------------------*/
#define MAX_P_FACTORS 6

typedef struct {
	uint32 p;
	uint32 r;
	uint8 log_p;
	uint8 max_roots;
	uint32 root_offsets[POLY_BATCH_SIZE + 1];
} sieve_prime_t;

typedef struct {
	sieve_prime_t *primes;
	uint32 num_primes;
	uint32 num_primes_alloc;
} sieve_prime_list_t;

typedef struct {
	uint8 *sieve_block;
	uint64 base;
	uint32 curr_offset;

	sieve_prime_list_t good_primes;
	sieve_prime_list_t bad_primes;

	uint32 num_small_roots;
	uint32 num_small_roots_alloc;
	uint32 *small_roots;

	uint32 max_roots;
	mpz_t *roots;

	mpz_t p, p2, nmodp2, m0, tmp1, tmp2;

	mpz_t accum[MAX_P_FACTORS + 1];

} sieve_fb_t;

void sieve_fb_init(sieve_fb_t *s, poly_batch_t *poly,
			uint32 factor_min, uint32 factor_max);

void sieve_fb_free(sieve_fb_t *s);

void sieve_fb_reset(sieve_fb_t *s, uint64 base);

uint64 sieve_fb_next(sieve_fb_t *s, poly_batch_t *poly, 
			p_batch_t *pbatch, uint64 limit,
			uint32 num_roots_min);

/*-----------------------------------------------------------------------*/

/* main search routines */

#define MAX_P_BITS 32
#define MAX_P ((uint64)1 << MAX_P_BITS)

void sieve_lattice(msieve_obj *obj, poly_batch_t *poly, 
			uint32 small_fb_max, uint32 large_fb_min, 
			uint32 large_fb_max, uint32 deadline);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_H_ */
