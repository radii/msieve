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

#ifndef _GNFS_POLY_POLY_SKEW_H_
#define _GNFS_POLY_POLY_SKEW_H_

#include "poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* external interface to skewed polynomial selector */

typedef void (*stage1_callback_t)(mpz_t high_coeff, mpz_t p, mpz_t m, 
				   double coeff_bound, void *extra);

typedef struct {
	mpz_t gmp_N;
	mpz_t gmp_high_coeff_begin;
	mpz_t gmp_high_coeff_end;
	uint32 degree;
	uint32 deadline;
	double norm_max;
	stage1_callback_t callback;
	void *callback_data;
} poly_stage1_t;

void poly_stage1_init(poly_stage1_t *data, 
			stage1_callback_t callback,
			void *callback_data);
void poly_stage1_free(poly_stage1_t *data);
uint32 poly_stage1_run(msieve_obj *obj, poly_stage1_t *data);


typedef void (*stage2_callback_t)(void *extra, uint32 deg,
				mpz_t * coeff1, mpz_t * coeff2,
				double skewness, double size_score,
				double root_score, double combined_score);

typedef struct {
	mpz_t gmp_N;
	uint32 degree;
	uint32 murphy_p_bound;
	double max_norm;
	double min_e;
	double min_e_bernstein;

	void *internal;

	stage2_callback_t callback;
	void *callback_data;
} poly_stage2_t;

void poly_stage2_init(poly_stage2_t *data,
		      stage2_callback_t callback,
		      void *callback_data);
void poly_stage2_free(poly_stage2_t *data);
void poly_stage2_run(poly_stage2_t *data, mpz_t a5, mpz_t p, 
			mpz_t d, double a3_bound);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_POLY_POLY_SKEW_H_ */
