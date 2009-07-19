#include "stage2.h"

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
static int
pol_expand(curr_poly_t *c, mpz_t gmp_N, double a3_bound)
{
	/* compute coefficients */
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[5]);
	mpz_sub(c->gmp_help3, gmp_N, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	if (mpz_cmp_ui(c->gmp_p, (mp_limb_t)1) == 0)
		mpz_set_ui(c->gmp_help2, (mp_limb_t)1);
	else {
		if (!mpz_invert(c->gmp_help2, c->gmp_d, c->gmp_p))
			return 0;
	}
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_a[4], c->gmp_help4, c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[4]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[3], c->gmp_help3, c->gmp_help4);
	mpz_fdiv_q(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
	mpz_mul(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[3], c->gmp_a[3], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[3]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[2], c->gmp_help3, c->gmp_help4);
	mpz_fdiv_q(c->gmp_a[2], c->gmp_a[2], c->gmp_p);
	mpz_mul(c->gmp_a[2], c->gmp_a[2], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[2]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_fdiv_q(c->gmp_a[1], c->gmp_help3, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[1], c->gmp_a[1], c->gmp_p);
	mpz_mul(c->gmp_a[1], c->gmp_a[1], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help3, c->gmp_help2);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[1], c->gmp_a[1], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_a[1]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_set(c->gmp_a[0], c->gmp_help3);

	mpz_fdiv_qr(c->gmp_help1, c->gmp_a[3], c->gmp_a[3], c->gmp_d);
	mpz_add(c->gmp_help2, c->gmp_a[3], c->gmp_a[3]);
	if (mpz_cmp(c->gmp_d, c->gmp_help2) < 0) {
		mpz_sub(c->gmp_a[3], c->gmp_a[3], c->gmp_d);
		mpz_add_ui(c->gmp_help1, c->gmp_help1, (mp_limb_t)1);
	}
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[4], c->gmp_a[4], c->gmp_help1);

	mpz_fdiv_qr(c->gmp_help1, c->gmp_a[2], c->gmp_a[2], c->gmp_d);
	mpz_add(c->gmp_help2, c->gmp_a[2], c->gmp_a[2]);
	if (mpz_cmp(c->gmp_d, c->gmp_help2) < 0) {
		mpz_sub(c->gmp_a[2], c->gmp_a[2], c->gmp_d);
		mpz_add_ui(c->gmp_help1, c->gmp_help1, (mp_limb_t)1);
	}
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[3], c->gmp_a[3], c->gmp_help1);

	mpz_set(c->gmp_lina[1], c->gmp_p);
	mpz_neg(c->gmp_lina[0], c->gmp_d);
#if 0
	uint32 i;
	gmp_printf("%+Zd\n", c->gmp_lina[0]);
	gmp_printf("%+Zd\n", c->gmp_lina[1]);
	for (i = 0; i <= 5; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);
#endif
	if (mpz_cmpabs_d(c->gmp_a[3], a3_bound) > 0)
		return 1;
	return 2;
}

/*-------------------------------------------------------------------------*/
static void
curr_poly_init(curr_poly_t *c)
{
	int i;

	mpz_init(c->gmp_p);
	mpz_init(c->gmp_d);
	for (i = 0; i < 2; i++)
		mpz_init(c->gmp_lina[i]);
	for (i = 0; i < 2; i++)
		mpz_init(c->gmp_linb[i]);
	for (i = 0; i < 6; i++)
		mpz_init(c->gmp_a[i]);
	for (i = 0; i < 6; i++)
		mpz_init(c->gmp_b[i]);
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
	for (i = 0; i < 2; i++)
		mpz_clear(c->gmp_lina[i]);
	for (i = 0; i < 2; i++)
		mpz_clear(c->gmp_linb[i]);
	for (i = 0; i < 6; i++)
		mpz_clear(c->gmp_a[i]);
	for (i = 0; i < 6; i++)
		mpz_clear(c->gmp_b[i]);
	mpz_clear(c->gmp_help1);
	mpz_clear(c->gmp_help2);
	mpz_clear(c->gmp_help3);
	mpz_clear(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_init(poly_stage2_t *data, 
		 stage2_callback_t callback,
		 void *callback_data)
{
	memset(data, 0, sizeof(poly_stage2_t));
	mpz_init(data->gmp_N);
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
poly_stage2_run(poly_stage2_t *data, mpz_t a5, mpz_t p, 
			mpz_t d, double a3_bound)
{
	double pol_norm;
	double alpha_proj;
	int status;
	stage2_curr_data_t *s;
	curr_poly_t *c;

	dd_precision_t precision = 0;
	uint32 precision_changed = 0;

	if (!dd_precision_is_ieee()) {
		precision_changed = 1;
		precision = dd_set_precision_ieee();
	}

	s = (stage2_curr_data_t *)(data->internal);
	if (s == NULL) {
		s = (stage2_curr_data_t *)xmalloc(sizeof(stage2_curr_data_t));
		curr_poly_init(&s->curr_poly);
		root_sieve_init(&s->root_sieve);
		assess_init(&s->assess);

		s->size_cutoff = data->min_e_bernstein / 
					pow(exp(-7.0), -2./6);
		data->internal = (void *)s;
	}

	c = &s->curr_poly;
	mpz_set(c->gmp_a[5], a5);
	mpz_set(c->gmp_p, p);
	mpz_set(c->gmp_d, d);
	status = pol_expand(c, data->gmp_N, a3_bound);
       	if (status != 2) {
		if (status == 0)
			fprintf(stderr, "expand failed\n");
		goto finished;
	}

	optimize_initial(data, &pol_norm);

	if (pol_norm <= s->size_cutoff)
		goto finished;

	stage2_root_score(5, c->gmp_a, 100, &alpha_proj, 1);
	root_sieve_run(data, alpha_proj);

finished:
	if (precision_changed)
		dd_clear_precision(precision);
}
