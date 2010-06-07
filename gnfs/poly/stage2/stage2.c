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

#include "stage2.h"

#if 0
#define CHECK
#endif

/*----------------------------------------------------------------------*/
void
assess_init(assess_t *a)
{
	integrate_init(&a->integ_aux, SIZE_EPS,
			double_exponential);
	dickman_init(&a->dickman_aux);
}

/*----------------------------------------------------------------------*/
void
assess_free(assess_t *a)
{
	integrate_free(&a->integ_aux);
	dickman_free(&a->dickman_aux);
}

/*-------------------------------------------------------------------------*/
#ifdef CHECK
static uint32
check_poly(curr_poly_t *c, mpz_t *coeffs, mpz_t lin0, 
		mpz_t gmp_N, uint32 degree) {

	uint32 i;

	mpz_set(c->gmp_help1, coeffs[degree]);
	mpz_set(c->gmp_help2, c->gmp_p);
	for (i = degree; i; i--) {
		mpz_mul(c->gmp_help1, c->gmp_help1, lin0);
		mpz_neg(c->gmp_help1, c->gmp_help1);
		mpz_addmul(c->gmp_help1, coeffs[i-1], c->gmp_help2);
		mpz_mul(c->gmp_help2, c->gmp_help2, c->gmp_p);
	}
	mpz_tdiv_r(c->gmp_help1, c->gmp_help1, gmp_N);
	if (mpz_cmp_ui(c->gmp_help1, (mp_limb_t)0) != 0) {
		printf("error: corrupt polynomial expand\n");
		return 0;
	}
	return 1;
}
#endif
/*-------------------------------------------------------------------------*/
static int
pol_expand(curr_poly_t *c, mpz_t gmp_N, mpz_t high_coeff,
		mpz_t gmp_p, mpz_t gmp_d, 
		double coeff_bound, uint32 degree)
{
	uint32 i;

	mpz_set(c->gmp_p, gmp_p);
	mpz_set(c->gmp_d, gmp_d);
	mpz_set(c->gmp_lina[1], gmp_p);
	mpz_neg(c->gmp_lina[0], gmp_d);

	if (mpz_cmp_ui(c->gmp_p, (mp_limb_t)1) == 0)
		mpz_set_ui(c->gmp_help1, (mp_limb_t)1);
	else {
		if (!mpz_invert(c->gmp_help1, gmp_d, gmp_p))
			return 0;
	}

	mpz_set(c->gmp_b[1], c->gmp_help1);
	for (i = 2; i < degree; i++)
		mpz_mul(c->gmp_b[i], c->gmp_b[i-1], c->gmp_help1);

	mpz_set(c->gmp_c[1], gmp_d);
	for (i = 2; i <= degree; i++)
		mpz_mul(c->gmp_c[i], c->gmp_c[i-1], gmp_d);

	mpz_set(c->gmp_a[degree], high_coeff);
	mpz_set(c->gmp_help2, gmp_N);

	for (i = degree - 1; (int32)i >= 0; i--) {

		mpz_mul(c->gmp_help3, c->gmp_a[i+1], c->gmp_c[i+1]);
		mpz_sub(c->gmp_help3, c->gmp_help2, c->gmp_help3);
		mpz_tdiv_q(c->gmp_help2, c->gmp_help3, gmp_p);

		if (i > 0) {
			mpz_tdiv_q(c->gmp_a[i], c->gmp_help2, c->gmp_c[i]);
			mpz_mul(c->gmp_help3, c->gmp_help2, c->gmp_b[i]);
			mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_a[i]);
			mpz_tdiv_r(c->gmp_help4, c->gmp_help3, gmp_p);

			if (mpz_sgn(c->gmp_help4) < 0)
				mpz_add(c->gmp_help4, c->gmp_help4, gmp_p);

			mpz_add(c->gmp_a[i], c->gmp_a[i], c->gmp_help4);
		}
	}
	mpz_set(c->gmp_a[0], c->gmp_help2);

	mpz_tdiv_q_2exp(c->gmp_help1, gmp_d, (mp_limb_t)1);
	for (i = 0; i < degree; i++) {
		while (mpz_cmpabs(c->gmp_a[i], c->gmp_help1) > 0) {
			if (mpz_sgn(c->gmp_a[i]) < 0) {
				mpz_add(c->gmp_a[i], c->gmp_a[i], gmp_d);
				mpz_sub(c->gmp_a[i+1], c->gmp_a[i+1], gmp_p);
			}
			else {
				mpz_sub(c->gmp_a[i], c->gmp_a[i], gmp_d);
				mpz_add(c->gmp_a[i+1], c->gmp_a[i+1], gmp_p);
			}
		}
	}

#if 0
	gmp_printf("%+Zd\n", c->gmp_lina[0]);
	gmp_printf("%+Zd\n", c->gmp_lina[1]);
	for (i = 0; i <= degree; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);

	printf("coeff ratio = %.5lf\n",
		fabs(mpz_get_d(c->gmp_a[degree-2])) / coeff_bound);
#endif

#ifdef CHECK
	if (check_poly(c, c->gmp_a, 
			c->gmp_lina[0], gmp_N, degree) != 1) {
		return 0;
	}
#endif

	if (mpz_cmpabs_d(c->gmp_a[degree - 2], coeff_bound) > 0) {
		return 1;
	}
	return 2;
}

/*-------------------------------------------------------------------------*/
static void
curr_poly_init(curr_poly_t *c)
{
	int i;

	mpz_init(c->gmp_p);
	mpz_init(c->gmp_d);
	for (i = 0; i < 2; i++) {
		mpz_init(c->gmp_lina[i]);
		mpz_init(c->gmp_linb[i]);
	}
	for (i = 0; i < MAX_POLY_DEGREE + 1; i++) {
		mpz_init(c->gmp_a[i]);
		mpz_init(c->gmp_b[i]);
		mpz_init(c->gmp_c[i]);
	}
	mpz_init(c->gmp_help1);
	mpz_init(c->gmp_help2);
	mpz_init(c->gmp_help3);
	mpz_init(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
static void
curr_poly_free(curr_poly_t *c)
{
	int i;

	mpz_clear(c->gmp_p);
	mpz_clear(c->gmp_d);
	for (i = 0; i < 2; i++) {
		mpz_clear(c->gmp_lina[i]);
		mpz_clear(c->gmp_linb[i]);
	}
	for (i = 0; i < MAX_POLY_DEGREE + 1; i++) {
		mpz_clear(c->gmp_a[i]);
		mpz_clear(c->gmp_b[i]);
		mpz_clear(c->gmp_c[i]);
	}
	mpz_clear(c->gmp_help1);
	mpz_clear(c->gmp_help2);
	mpz_clear(c->gmp_help3);
	mpz_clear(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_init(poly_stage2_t *data, 
		 msieve_obj *obj,
		 stage2_callback_t callback,
		 void *callback_data)
{
	memset(data, 0, sizeof(poly_stage2_t));
	mpz_init(data->gmp_N);
	data->obj = obj;
	data->murphy_p_bound = 2000;
	data->callback = callback;
	data->callback_data = callback_data;
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_free(poly_stage2_t *data)
{
	mpz_clear(data->gmp_N);

	if (data->internal) {
		stage2_curr_data_t *s = (stage2_curr_data_t *)(data->internal);
		curr_poly_free(&s->curr_poly);
		root_sieve_free(&s->root_sieve);
		assess_free(&s->assess);

		free(data->internal);
		data->internal = NULL;
	}
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_run(poly_stage2_t *data, mpz_t high_coeff, mpz_t p, 
			mpz_t d, double coeff_bound)
{
	double pol_norm;
	double alpha_proj;
	int status;
	stage2_curr_data_t *s;
	curr_poly_t *c;
	uint32 degree;

	dd_precision_t precision = 0;
	uint32 precision_changed = 0;

	if (!dd_precision_is_ieee()) {
		precision_changed = 1;
		precision = dd_set_precision_ieee();
	}

	degree = data->degree;
	s = (stage2_curr_data_t *)(data->internal);
	if (s == NULL) {
		s = (stage2_curr_data_t *)xmalloc(sizeof(stage2_curr_data_t));
		curr_poly_init(&s->curr_poly);
		root_sieve_init(&s->root_sieve);
		assess_init(&s->assess);
		data->internal = (void *)s;
	}

	c = &s->curr_poly;
	status = pol_expand(c, data->gmp_N, high_coeff,
				p, d, coeff_bound, degree);
       	if (status != 2) {
		if (status == 0)
			fprintf(stderr, "expand failed\n");
		goto finished;
	}

	optimize_initial(data, &pol_norm);

#ifdef CHECK
	if (check_poly(c, c->gmp_a, c->gmp_lina[0],
			data->gmp_N, degree) != 1)
		goto finished;

	printf("%le %le\n", pol_norm, data->max_norm);
#endif
	if (pol_norm > data->max_norm)
		goto finished;

	stage2_root_score(degree, c->gmp_a, 100, &alpha_proj, 1);
	root_sieve_run(data, pol_norm, alpha_proj);

finished:
	if (precision_changed)
		dd_clear_precision(precision);
}
