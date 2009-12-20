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
#include <cuda_xface.h>

#ifdef __cplusplus
extern "C" {
#endif

#define POLY_BATCH_SIZE 40

#define MAX_POLYSELECT_DEGREE 6

#if MAX_POLY_DEGREE < MAX_POLYSELECT_DEGREE
#error "supported poly degree must be at least 6"
#endif

#define MULTIPLIER 60	/* 2*2*3*5 */

/* 96-bit integers */

typedef struct {
	uint32 w[3];
} uint96;

/* 128-bit integers */

typedef struct {
	uint32 w[4];
} uint128;

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
	mpz_t high_coeff; 
	mpz_t trans_N;
	mpz_t trans_m0;

	double coeff_max;
	double p_size_max;

	double sieve_size;
	mpz_t mp_sieve_size;
} curr_poly_t;

typedef struct {

	uint32 degree;
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

#define MAX_P_FACTORS 7
#define MAX_ROOTS 128

#define P_SEARCH_DONE ((uint64)(-1))

typedef struct {
	uint32 p;
	uint32 r;
	float fp_log_p;
	uint8 log_p;
	uint8 num_roots[POLY_BATCH_SIZE];
	uint32 roots[POLY_BATCH_SIZE][MAX_POLYSELECT_DEGREE];
} sieve_prime_t;

typedef struct {
	sieve_prime_t *primes;
	uint32 num_primes;
	uint32 num_primes_alloc;
} sieve_prime_list_t;

typedef struct {
	uint32 num_factors;
	float log_prod;
	uint64 prod;
} ss_t;

typedef struct {
	ss_t *list;
	uint32 num_entries;
	uint32 num_entries_alloc;
} subset_sum_t;

typedef struct {
	subset_sum_t curr_list;
	subset_sum_t new_list;

	uint32 curr_entry;
	uint32 next_prime;
	uint32 num_primes;
	float log_p_min;
	float log_p_max;
} p_enum_t;

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

	p_enum_t p_enum;

	mpz_t p, p2, m0, nmodp2, tmp1, tmp2;
	mpz_t accum[MAX_P_FACTORS + 1];
	mpz_t roots[MAX_ROOTS];
} sieve_fb_t;

void sieve_fb_init(sieve_fb_t *s, poly_search_t *poly,
			uint32 factor_min, uint32 factor_max,
			uint32 fb_roots_min, uint32 fb_roots_max);

void sieve_fb_free(sieve_fb_t *s);

void sieve_fb_reset(sieve_fb_t *s, uint64 p_min, uint64 p_max,
			uint32 num_roots_min, uint32 num_roots_max);

typedef void (*root_callback)(uint64 p, uint32 num_roots, 
				uint32 which_poly, mpz_t *roots, 
				void *extra);

uint64 sieve_fb_next(sieve_fb_t *s, 
			poly_search_t *poly, 
			root_callback callback,
			void *extra);

/*-----------------------------------------------------------------------*/

typedef struct {
	uint32 fill_p;
	void *p_array;
	void *q_array;

	CUdeviceptr gpu_p_array;
	CUdeviceptr gpu_q_array;
	CUdeviceptr gpu_found_array;
	void *found_array;
	uint32 found_array_size;
	void *p_marshall;
	void *q_marshall;

	poly_search_t *poly;

	time_t start_time;
	uint32 deadline;
	uint32 num_tests;
	uint32 tests_per_block;
} lattice_fb_t;

/* lower-level sieve routines */

uint32
sieve_lattice_gpu_deg46_64(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

uint32
sieve_lattice_gpu_deg5_64(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

uint32
sieve_lattice_gpu_deg5_96(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

uint32
sieve_lattice_gpu_deg5_128(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

uint32
sieve_lattice_gpu_deg6_96(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

uint32
sieve_lattice_gpu_deg6_128(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint64 small_p_min, uint64 small_p_max, 
		uint64 large_p_min, uint64 large_p_max,
		gpu_info_t *gpu_info, CUfunction gpu_kernel);

void
handle_collision(poly_search_t *poly, uint32 which_poly,
		uint64 p, uint128 proot, uint128 res, uint64 q);

/* main search routine */

void sieve_lattice(msieve_obj *obj, poly_search_t *poly, 
			uint32 small_fb_max, uint32 large_fb_min, 
			uint32 large_fb_max, gpu_info_t *gpu_info,
			CUmodule gpu_module64, CUmodule gpu_module96, 
			CUmodule gpu_module128, uint32 deadline);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_H_ */
