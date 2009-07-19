#ifndef _STAGE2_H_
#define _STAGE2_H_

#include <poly_skew.h>

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------*/
/* data used in the current polynomial */

typedef struct {
	mpz_t gmp_a[6];
	mpz_t gmp_b[6];
	mpz_t gmp_lina[2];
	mpz_t gmp_linb[2];
	mpz_t gmp_help1;
	mpz_t gmp_help2;
	mpz_t gmp_help3;
	mpz_t gmp_help4;
	mpz_t gmp_p;
	mpz_t gmp_d;
} curr_poly_t;

/*-----------------------------------------------------------------------*/
/* data for rating polynomial yield */

typedef struct {
	integrate_t integ_aux;
	dickman_t dickman_aux;
} assess_t;

void assess_init(assess_t *a);

void assess_free(assess_t *a);

uint32 stage2_root_score(uint32 deg1, mpz_t *coeff1, 
       			uint32 prime_bound, double *score,
       			uint32 projective_only);

/*-----------------------------------------------------------------------*/
/* routines for optimizing polynomials */

void optimize_initial(poly_stage2_t *data, double *pol_norm);

void optimize_final(int64 x, int y, poly_stage2_t *data);

double optimize_basic(dpoly_t *apoly, double *best_skewness,
				double *best_translation);

/*-----------------------------------------------------------------------*/
/* data for the root sieve */

typedef struct {
	uint16 resclass;
	uint16 start;
	uint16 step;
} sieve_root_t;

typedef struct {
	uint32 power;
	uint32 num_roots;
	double root_contrib;
	uint16 sieve_contrib;
	sieve_root_t *roots;
} sieve_power_t;

typedef struct {
	uint32 prime;
	uint32 num_powers;
	sieve_power_t *powers;

	uint32 contrib_array_size;
	uint32 contrib_array_offset;
	uint16 *contrib_array;
} sieve_prime_t;

typedef struct {
	int64 x;
	int32 y;
	uint16 score;
} rotation_t;

typedef struct {
	uint32 num_entries;
	uint32 max_entries;
	rotation_t *entries;

	rotation_t cutoffs[2];
	uint32 default_cutoff;
	void *extra;
} root_heap_t;

typedef struct {
	uint32 num_primes;
	sieve_prime_t *primes;

	double sieve_bias; 
	double random_root_score;

	uint16 *sieve_block;

	root_heap_t root_heap;
	root_heap_t lattice_heap;
	root_heap_t tmp_lattice_heap;
} root_sieve_t;

void root_sieve_init(root_sieve_t *rs);
void root_sieve_free(root_sieve_t *rs);
void root_sieve_run(poly_stage2_t *data, double alpha_proj);

/*-------------------------------------------------------------------------*/

/* data for optimizing a single (a5, p, d) triplet */

typedef struct {
	curr_poly_t curr_poly;
	root_sieve_t root_sieve;
	assess_t assess;
	double size_cutoff;
} stage2_curr_data_t;


#ifdef __cplusplus
}
#endif

#endif /* !_STAGE2_H_ */
