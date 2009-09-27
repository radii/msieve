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

#include "stage1.h"

/*------------------------------------------------------------------------*/
void 
stage1_bounds_init(bounds_t *bounds, poly_stage1_t *data)
{
	mpz_init_set(bounds->gmp_high_coeff_begin, 
			data->gmp_high_coeff_begin);
	mpz_init_set(bounds->gmp_high_coeff_end, 
			data->gmp_high_coeff_end);
	bounds->norm_max = data->norm_max;
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_free(bounds_t *bounds)
{
	mpz_clear(bounds->gmp_high_coeff_begin);
	mpz_clear(bounds->gmp_high_coeff_end);
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_update(bounds_t *bounds, double N, 
			double high_coeff, uint32 degree)
{
	double skewness_min = 1.0;

	switch (degree) {
	case 4:
		skewness_min = sqrt(pow(N / high_coeff, 1./4.) / 
					bounds->norm_max);
		bounds->coeff_max = bounds->norm_max;
		break;

	case 5:
		skewness_min = pow(pow(N / high_coeff, 1./5.) / 
					bounds->norm_max, 2./3.);
		bounds->coeff_max = bounds->norm_max / sqrt(skewness_min);
		break;

	case 6:
		skewness_min = sqrt(pow(N / high_coeff, 1./6.) / 
					bounds->norm_max);
		bounds->coeff_max = bounds->norm_max / skewness_min;
		break;

	default:
		printf("error: unexpected poly degree %d\n", degree);
		exit(-1);
	}

	bounds->p_size_max = bounds->coeff_max / skewness_min;
}

/*------------------------------------------------------------------------*/
void
poly_search_init(poly_search_t *poly, poly_stage1_t *data)
{
	mpz_init(poly->high_coeff);
	mpz_init(poly->trans_N);
	mpz_init(poly->trans_m0);
	mpz_init_set(poly->N, data->gmp_N);
	mpz_init(poly->m0);
	mpz_init(poly->p);
	mpz_init(poly->tmp1);
	mpz_init(poly->tmp2);
	mpz_init(poly->tmp3);
	mpz_init(poly->tmp4);
	mpz_init(poly->tmp5);

	poly->degree = data->degree;
	poly->callback = data->callback;
	poly->callback_data = data->callback_data;
}

/*------------------------------------------------------------------------*/
void
poly_search_free(poly_search_t *poly)
{
	mpz_clear(poly->high_coeff);
	mpz_clear(poly->trans_N);
	mpz_clear(poly->trans_m0);
	mpz_clear(poly->N);
	mpz_clear(poly->m0);
	mpz_clear(poly->p);
	mpz_clear(poly->tmp1);
	mpz_clear(poly->tmp2);
	mpz_clear(poly->tmp3);
	mpz_clear(poly->tmp4);
	mpz_clear(poly->tmp5);
}

/*------------------------------------------------------------------------*/
static void
search_coeffs_core(msieve_obj *obj, poly_search_t *poly, uint32 deadline)
{
	uint32 i;
	uint32 degree = poly->degree;
	uint32 mult = 0;

	switch (degree) {
	case 4: mult = 4 * 4 * 4 * 4; break;
	case 5: mult = 5 * 5 * 5 * 5 * 5; break;
	case 6: mult = 6 * 6 * 6 * 6 * 6 * 6; break;
	}

	mpz_mul_ui(poly->trans_N, poly->N, (mp_limb_t)mult);
	for (i = 0; i < degree - 1; i++)
		mpz_mul(poly->trans_N, poly->trans_N, poly->high_coeff);

	mpz_root(poly->trans_m0, poly->trans_N, (mp_limb_t)degree);
	mpz_tdiv_q(poly->m0, poly->N, poly->high_coeff);
	mpz_root(poly->m0, poly->m0, (mp_limb_t)degree);

	poly->sieve_size = poly->coeff_max / mpz_get_d(poly->m0) * 
			poly->p_size_max * poly->p_size_max / degree;

	sieve_lattice(obj, poly, 2000, 2001, 100000, deadline);
}

/*------------------------------------------------------------------------*/
static void
search_coeffs(msieve_obj *obj, poly_search_t *poly, 
		bounds_t *bounds, uint32 deadline)
{
	mpz_t curr_high_coeff;
	double dn = mpz_get_d(poly->N);
	uint32 digits = mpz_sizeinbase(poly->N, 10);
	double start_time = get_cpu_time();
	uint32 deadline_per_coeff = 800;

	if (digits <= 100)
		deadline_per_coeff = 5;
	else if (digits <= 105)
		deadline_per_coeff = 20;
	else if (digits <= 110)
		deadline_per_coeff = 30;
	else if (digits <= 120)
		deadline_per_coeff = 50;
	else if (digits <= 130)
		deadline_per_coeff = 100;
	else if (digits <= 140)
		deadline_per_coeff = 200;
	else if (digits <= 150)
		deadline_per_coeff = 400;
	printf("deadline: %u seconds per coefficient\n", deadline_per_coeff);

	mpz_init(curr_high_coeff);
	mpz_fdiv_q_ui(curr_high_coeff, 
			bounds->gmp_high_coeff_begin, 
			(mp_limb_t)MULTIPLIER);
	mpz_mul_ui(curr_high_coeff, curr_high_coeff, 
			(mp_limb_t)MULTIPLIER);
	if (mpz_cmp(curr_high_coeff, 
			bounds->gmp_high_coeff_begin) < 0) {
		mpz_add_ui(curr_high_coeff, curr_high_coeff, 
				(mp_limb_t)MULTIPLIER);
	}

	while (1) {
		if (mpz_cmp(curr_high_coeff, bounds->gmp_high_coeff_end) > 0)
			break;

		stage1_bounds_update(bounds, dn, 
					mpz_get_d(curr_high_coeff),
					poly->degree);

		mpz_set(poly->high_coeff, curr_high_coeff);
		poly->p_size_max = bounds->p_size_max;
		poly->coeff_max = bounds->coeff_max;

		search_coeffs_core(obj, poly, deadline_per_coeff);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

		if (deadline) {
			double curr_time = get_cpu_time();
			double elapsed = curr_time - start_time;

			if (elapsed > deadline)
				break;
		}

		mpz_add_ui(curr_high_coeff, curr_high_coeff, 
				(mp_limb_t)MULTIPLIER);
	}

	mpz_clear(curr_high_coeff);
}

/*------------------------------------------------------------------------*/
void
poly_stage1_init(poly_stage1_t *data,
		 stage1_callback_t callback, void *callback_data)
{
	memset(data, 0, sizeof(poly_stage1_t));
	mpz_init_set_ui(data->gmp_N, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_begin, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_end, (mp_limb_t)0);
	data->callback = callback;
	data->callback_data = callback_data;
}

/*------------------------------------------------------------------------*/
void
poly_stage1_free(poly_stage1_t *data)
{
	mpz_clear(data->gmp_N);
	mpz_clear(data->gmp_high_coeff_begin);
	mpz_clear(data->gmp_high_coeff_end);
}

/*------------------------------------------------------------------------*/
uint32
poly_stage1_run(msieve_obj *obj, poly_stage1_t *data)
{
	bounds_t bounds;
	poly_search_t poly;

	stage1_bounds_init(&bounds, data);
	poly_search_init(&poly, data);
	search_coeffs(obj, &poly, &bounds, data->deadline);

	poly_search_free(&poly);
	stage1_bounds_free(&bounds);
	return 1;
}
