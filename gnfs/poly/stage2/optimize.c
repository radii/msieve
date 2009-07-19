#include "stage2.h"

/*-------------------------------------------------------------------------*/
static double
ifs(double *a, uint32 degree, double s)
{
	double a0, a1, a2, a3, a4, a5, a6;
	double s2, s3, s4, s5, s6;
	double norm;

	if (s < 1)
		return 1e100;

	a0 = a[0];
	a1 = a[1] * s;
	s2 = s * s;
	a2 = a[2] * s2;
	s3 = s2 * s;
	a3 = a[3] * s3;
	s4 = s3 * s;
	a4 = a[4] * s4;

	switch (degree) {
	case 4:
		norm = 1.0 / 9.0 * (a4 * a4 + a0 * a0) + 
		       2.0 / 21.0 * (a4 * a2 + a2 * a0) + 
		       1.0 / 21.0 * (a3 * a3 + a1 * a1) + 
		       2.0 / 25.0 * (a4 * a0 + a3 * a1) + 
		       1.0 / 25.0 * a2 * a2;
		return norm / s4;

	case 5:
		s5 = s4 * s;
		a5 = a[5] * s5;
		norm = 1.0 / 11.0 * (a5 * a5 + a0 * a0) + 
		       2.0 / 27.0 * (a5 * a3 + a2 * a0) + 
		       1.0 / 27.0 * (a4 * a4 + a1 * a1) +
		       1.0 / 35.0 * (a3 * a3 + a2 * a2) +
		       2.0 / 35.0 * (a5 * a1 + a4 * a2 + a4 * a0 + a3 * a1);
		return norm / s5;

	case 6:
		s5 = s4 * s;
		a5 = a[5] * s5;
		s6 = s5 * s;
		a6 = a[6] * s6;
		norm = 1.0 / 13.0 * (a6 * a6 + a0 * a0) + 
		       2.0 / 33.0 * (a6 * a4 + a2 * a0) + 
		       1.0 / 33.0 * (a5 * a5 + a1 * a1) +
		       1.0 / 45.0 * (a4 * a4 + a2 * a2) + 
		       2.0 / 45.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1) + 
		       2.0 / 49.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 
		       1.0 / 49.0 * a3 * a3;
		return norm / s6;
	}

	return 1e100;
}

/*----------------------------------------------------------------------*/
uint32
stage2_root_score(uint32 deg1, mpz_t *coeff1, 
		 uint32 prime_bound, double *score,
		 uint32 projective_only)
{
	uint32 i;
	mp_poly_t apoly;

	memset(&apoly, 0, sizeof(mp_poly_t));

	apoly.degree = deg1;
	for (i = 0; i <= deg1; i++) {
		gmp2mp(coeff1[i], &apoly.coeff[i].num);
		apoly.coeff[i].sign = POSITIVE;
		if (mpz_sgn(coeff1[i]) < 0)
			apoly.coeff[i].sign = NEGATIVE;
	}

	if (projective_only)
		return analyze_poly_roots_projective(&apoly, 
						prime_bound, score);
	else
		return analyze_poly_roots(&apoly, prime_bound, score);
}

static const double xlate_weights[6+1][6+1] = {
	{  1.,  1.,  1.,  1.,  1.,  1.,  1.},
	{  0.,  1.,  2.,  3.,  4.,  5.,  6.},
	{  0.,  0.,  1.,  3.,  6., 10., 15.},
	{  0.,  0.,  0.,  1.,  4., 10., 20.},
	{  0.,  0.,  0.,  0.,  1.,  5., 15.},
	{  0.,  0.,  0.,  0.,  0.,  1.,  6.},
	{  0.,  0.,  0.,  0.,  0.,  0.,  1.},
};

/*-------------------------------------------------------------------------*/
static void
translate_d(double *c, double *d, uint32 deg, double x)
{
	uint32 i, j;

	for (i = 0, x = -x; i < deg; i++) {

		const double *weights = xlate_weights[i];
		double accum = d[deg] * weights[deg];

		for (j = deg; j > i; j--) {
			accum = accum * x + d[j-1] * weights[j-1];
		}
		c[i] = accum;
	}
	c[i] = d[i];
}

/*-------------------------------------------------------------------------*/
static void
translate_dd(dd_t *c, dd_t *d, uint32 deg, double x)
{
	uint32 i, j;

	for (i = 0, x = -x; i < deg; i++) {

		const double *weights = xlate_weights[i];
		dd_t accum = dd_mul_d(d[deg], weights[deg]);

		for (j = deg; j > i; j--) {
			accum = dd_add_dd(dd_mul_d(accum, x),
					  dd_mul_d(d[j-1], weights[j-1]));
		}
		c[i] = accum;
	}
	c[i] = d[i];
}

/*-------------------------------------------------------------------------*/
static void
translate_gmp(curr_poly_t *c, mpz_t *gmp_c, uint32 deg,
		mpz_t *lin, int64 k)
{
	/* note that unlike the other translation functions, 
	   gmp_c[] and lin[] are overwritten */

	uint32 i, j;

	int64_2gmp(-k, c->gmp_help1);

	for (i = 0; i < deg; i++) {

		const double *weights = xlate_weights[i];

		mpz_set_d(c->gmp_help2, weights[deg]);
		mpz_mul(c->gmp_help2, c->gmp_help2, gmp_c[deg]);

		for (j = deg; j > i; j--) {
			mpz_set_d(c->gmp_help3, weights[j-1]);
			mpz_mul(c->gmp_help3, c->gmp_help3, gmp_c[j-1]);
			mpz_addmul(c->gmp_help3, c->gmp_help2, c->gmp_help1);
			mpz_swap(c->gmp_help2, c->gmp_help3);
		}
		mpz_swap(gmp_c[i], c->gmp_help2);
	}

	mpz_addmul(lin[0], lin[1], c->gmp_help1);
}

/*-------------------------------------------------------------------------*/
#define TRANSLATE_SIZE 0
#define SKEWNESS 1
#define ROTATE0 2
#define ROTATE1 3
#define ROTATE2 4

typedef struct {
	uint32 rotate_dim;
	dpoly_t *drpoly;
	dpoly_t *dapoly;

	ddpoly_t *rpoly;
	ddpoly_t *apoly;
	integrate_t *integ_aux;

	dickman_t *dickman_aux;
	double root_score_r;
	double root_score_a;
} opt_data_t;

static double poly_xlate_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	dpoly_t *apoly = opt->dapoly;
	double translated[MAX_POLY_DEGREE + 1];
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double s = v[SKEWNESS];

	translate_d(translated, apoly->coeff, apoly->degree, t);

	return ifs(translated, apoly->degree, s);
}

static double poly_rotate_callback(double *v, void *extra)
{
	uint32 i;
	opt_data_t *opt = (opt_data_t *)extra;
	dpoly_t apoly = *(opt->dapoly);
	double translated[MAX_POLY_DEGREE + 1];
	double s = v[SKEWNESS];
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double r0 = opt->drpoly->coeff[0];
	double r1 = opt->drpoly->coeff[1];

	if (s < 1.0)
		return 1e100;

	for (i = 0; i <= opt->rotate_dim; i++) {
		double c = floor(v[ROTATE0 + i] + 0.5);
		apoly.coeff[i] += r0 * c;
		apoly.coeff[i+1] += r1 * c;
	}

	translate_d(translated, apoly.coeff, apoly.degree, t);
	return ifs(translated, apoly.degree, s);
}

static double poly_murphy_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	ddpoly_t *apoly = opt->apoly;
	ddpoly_t *rpoly = opt->rpoly;
	ddpoly_t new_rpoly, new_apoly;
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double s = v[SKEWNESS];
	double score = 0;

	if (s < 0.5)
		return 1e50;

	new_rpoly.degree = 1;
	new_rpoly.coeff[1] = rpoly->coeff[1];
	new_rpoly.coeff[0] = dd_sub_dd(rpoly->coeff[0],
					dd_mul_d(rpoly->coeff[1], t));

	new_apoly.degree = apoly->degree;
	translate_dd(new_apoly.coeff, apoly->coeff, apoly->degree, t);

	analyze_poly_murphy(opt->integ_aux, opt->dickman_aux,
				&new_rpoly, opt->root_score_r,
				&new_apoly, opt->root_score_a,
				s, &score);

	return -score;
}

/*-------------------------------------------------------------------------*/
void
optimize_initial(poly_stage2_t *data, double *pol_norm)
{
	stage2_curr_data_t *curr_data = 
			(stage2_curr_data_t *)(data->internal);
	curr_poly_t *c = &curr_data->curr_poly;
	uint32 deg = data->degree;
	uint32 rotate_dim = deg - 4;
	opt_data_t opt_data;
	uint32 i;
	double best[MAX_VARS];
	double score, last_score;
	double s0, s1;
	dpoly_t rpoly, apoly;

	opt_data.rotate_dim = rotate_dim;
	opt_data.drpoly = &rpoly;
	opt_data.dapoly = &apoly;

	rpoly.degree = 1;
	for (i = 0; i <= 1; i++) {
		rpoly.coeff[i] = mpz_get_d(c->gmp_lina[i]);
	}
	apoly.degree = deg;
	for (i = 0; i <= deg; i++) {
		apoly.coeff[i] = mpz_get_d(c->gmp_a[i]);
	}

	s0 = s1 = 1.0;
	if (apoly.coeff[deg-1] != 0)
		s0 = fabs(apoly.coeff[deg-1] / 
			  apoly.coeff[deg]);
	if (apoly.coeff[deg-2] != 0)
		s1 = sqrt(fabs(apoly.coeff[deg-2] / 
			       apoly.coeff[deg]));
	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = MAX(s0, s1);
	best[ROTATE0] = 0;
	best[ROTATE1] = 0;
	best[ROTATE2] = 0;

	score = 1e100;
	do {
		last_score = score;
		score = minimize(best, rotate_dim + 3, 1e-5, 40, 
				poly_rotate_callback, &opt_data);

		for (i = 0; i <= rotate_dim; i++) {
			double ci = floor(best[ROTATE0 + i] + 0.5);
			mpz_set_d(c->gmp_help1, ci);
			mpz_addmul(c->gmp_a[i+1], c->gmp_help1, c->gmp_p);
			mpz_submul(c->gmp_a[i], c->gmp_help1, c->gmp_d);
		}
		translate_gmp(c, c->gmp_a, deg, c->gmp_lina,
				(int64)(best[TRANSLATE_SIZE] + 0.5));

		mpz_neg(c->gmp_d, c->gmp_lina[0]);
		for (i = 0; i <= 1; i++)
			rpoly.coeff[i] = mpz_get_d(c->gmp_lina[i]);
		for (i = 0; i <= deg; i++)
			apoly.coeff[i] = mpz_get_d(c->gmp_a[i]);
		best[TRANSLATE_SIZE] = 0;
		best[ROTATE0] = 0;
		best[ROTATE1] = 0;
		best[ROTATE2] = 0;

	} while (fabs(score - last_score) > .001 * fabs(score));

	*pol_norm = sqrt(fabs(score));
#if 0
	printf("norm %.7e skew %lf\n", *pol_norm, best[SKEWNESS]);
	for (i = 0; i < 2; i++)
		gmp_printf("%+Zd\n", c->gmp_lina[i]);
	for (i = 0; i <= deg; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);
#endif
}

/*-------------------------------------------------------------------------*/
double
optimize_basic(dpoly_t *apoly, double *best_skewness,
		double *best_translation)
{
	uint32 d = apoly->degree;
	opt_data_t opt_data;
	double best[MAX_VARS];
	double score;
	double s0 = 1.0;
	double s1 = 1.0;

	opt_data.dapoly = apoly;
	if (apoly->coeff[d-1] != 0)
		s0 = fabs(apoly->coeff[d-1] / apoly->coeff[d]);
	if (apoly->coeff[d-2] != 0)
		s1 = sqrt(fabs(apoly->coeff[d-2] / apoly->coeff[d]));

	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = MAX(s0, s1);

	score = minimize(best, 2, 1e-5, 40, poly_xlate_callback, &opt_data);

	*best_translation = floor(best[TRANSLATE_SIZE] + 0.5);
	*best_skewness  = best[SKEWNESS];
	return sqrt(score);
}

/*-------------------------------------------------------------------------*/
static void
optimize_final_core(curr_poly_t *c, assess_t *assess, uint32 deg,
			double root_score, double *best_score_out,
			double *best_skewness_out)
{
	uint32 i;
	signed_mp_t tmp;
	opt_data_t opt_data;
	double best[MAX_VARS];
	double score;
	ddpoly_t rpoly, apoly;

	opt_data.rpoly = &rpoly;
	opt_data.apoly = &apoly;
	opt_data.integ_aux = &assess->integ_aux;
	opt_data.dickman_aux = &assess->dickman_aux;
	opt_data.root_score_r = 0.0;
	opt_data.root_score_a = root_score;
	rpoly.degree = 1;
	apoly.degree = deg;

	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = sqrt(fabs(mpz_get_d(c->gmp_b[deg-2]) / 
				mpz_get_d(c->gmp_b[deg])));
	best[SKEWNESS] = MAX(1.0, best[SKEWNESS]);

	for (i = 0; i <= 1; i++) {
		gmp2mp(c->gmp_linb[i], &tmp.num);
		tmp.sign = mpz_sgn(c->gmp_linb[i]) >= 0 ? 
					POSITIVE : NEGATIVE;
		rpoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}
	for (i = 0; i <= deg; i++) {
		gmp2mp(c->gmp_b[i], &tmp.num);
		tmp.sign = mpz_sgn(c->gmp_b[i]) >= 0 ? 
					POSITIVE : NEGATIVE;
		apoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}

	score = minimize(best, 2, 1e-5, 40, 
			poly_murphy_callback, &opt_data);

	translate_gmp(c, c->gmp_b, deg, c->gmp_linb, 
			(int64)(best[TRANSLATE_SIZE] + 0.5));

	*best_score_out = fabs(score);
	*best_skewness_out = best[SKEWNESS];
}

/*-------------------------------------------------------------------------*/
static void
get_bernstein_score(curr_poly_t *c, assess_t *assess,
		uint32 deg, double root_score, double *eptr)
{
	uint32 i;
	double size_score;
	ddpoly_t rpoly, apoly;
	signed_mp_t tmp;

	root_score = pow(exp(root_score), -2./(deg + 1));

	rpoly.degree = 1;
	for (i = 0; i <= 1; i++) {
		gmp2mp(c->gmp_linb[i], &tmp.num);
		tmp.sign = mpz_sgn(c->gmp_linb[i]) >= 0 ? POSITIVE : NEGATIVE;
		rpoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}

	apoly.degree = deg;
	for (i = 0; i <= deg; i++) {
		gmp2mp(c->gmp_b[i], &tmp.num);
		tmp.sign = mpz_sgn(c->gmp_b[i]) >= 0 ? POSITIVE : NEGATIVE;
		apoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}

	analyze_poly_size(&assess->integ_aux, 
			&rpoly, &apoly, &size_score);
	if (size_score == 0.0)
		printf("error: size score computation failed\n");

	*eptr = size_score * root_score;
}

/*-------------------------------------------------------------------------*/
void
optimize_final(int64 x, int y, poly_stage2_t *data)
{
	uint32 i;
	uint32 deg = data->degree;
	double alpha, skewness, bscore, combined_score;
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	assess_t *assess = &s->assess;

	for (i = 0; i <= 1; i++)
		mpz_set(c->gmp_linb[i], c->gmp_lina[i]);

	for (i = 0; i <= deg; i++)
		mpz_set(c->gmp_b[i], c->gmp_a[i]);

	mpz_set_si(c->gmp_help1, y);
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_b[2], c->gmp_b[2], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_b[1], c->gmp_b[1], c->gmp_help1);

	int64_2gmp(x, c->gmp_help1);
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_b[1], c->gmp_b[1], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_b[0], c->gmp_b[0], c->gmp_help1);

	if (stage2_root_score(deg, c->gmp_b, data->murphy_p_bound, &alpha, 0))
		return;

	if (alpha > -4.5)
		return;

	get_bernstein_score(c, assess, deg, alpha, &bscore);

#if 0
	printf("%.0lf %d %lf %le\n", (double)x, y, alpha, bscore); 
	fflush(stdout);
#endif

	if (bscore > data->min_e_bernstein) {

		optimize_final_core(c, assess, deg, alpha, 
				&combined_score, &skewness);

#if 0
		printf("combined %le ratio %lf\n", combined_score,
					combined_score / bscore);
#endif
		if (combined_score > data->min_e) {
			data->callback(data->callback_data, deg, c->gmp_b, 
					c->gmp_linb, skewness, bscore,
					alpha, combined_score);
		}
	}
}
