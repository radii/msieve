#include "stage2.h"

/*-------------------------------------------------------------------------*/
static double
ifs(double *coeff, double skewness)
{				/* degree=5 */
	double sq[6], d, s, res;
	int i, j, k;

	if (skewness < 1)
		return 1e100;

	for (i = 0; i < 6; i++)
		sq[i] = coeff[i] * coeff[i];
	for (i = 0; i < 4; i++) {
		d = 2 * coeff[i];
		k = i + 1;
		for (j = i + 2; j < 6; j += 2, k++)
			sq[k] += d * coeff[j];
	}
	s = skewness;
	res = 0.;
	res += s * sq[3] / 35.;
	res += sq[2] / s / 35.;
	s *= (skewness * skewness);
	res += s * sq[4] / 27.;
	res += sq[1] / s / 27.;
	s *= (skewness * skewness);
	res += s * sq[5] / 11.;
	res += sq[0] / s / 11.;
	return res;
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

/*-------------------------------------------------------------------------*/
static void
translate_dbl(double *c, dd_t *d, double x)
{
	double d0 = d[0].hi;
	double d1 = d[1].hi;
	double d2 = d[2].hi;
	double d3 = d[3].hi;
	double d4 = d[4].hi;
	double d5 = d[5].hi;

	c[5] = d5;
	c[4] = d4 + x * (-5 * d5);
	c[3] = d3 + x * (-4 * d4 + x * 10 * d5);
	c[2] = d2 + x * (-3 * d3 + x * (6 * d4 - x * 10 * d5));
	c[1] = d1 + x * (-2 * d2 + x * (3 * d3 + x * (-4 * d4 + x * 5 * d5)));
	c[0] = d0 + x * (-d1 + x * (d2 + x * (-d3 + x * (d4 - x * d5))));
}

/*-------------------------------------------------------------------------*/
static void
translate_dd(dd_t *c, dd_t *d, double x)
{
	dd_t d0 = d[0];
	dd_t d1 = d[1];
	dd_t d2 = d[2];
	dd_t d3 = d[3];
	dd_t d4 = d[4];
	dd_t d5 = d[5];
	dd_t t;

	t = dd_sub_dd(d4, dd_mul_d(d5, x));
	t = dd_sub_dd(dd_mul_d(t, x), d3);
	t = dd_add_dd(d2, dd_mul_d(t, x));
	t = dd_sub_dd(dd_mul_d(t, x), d1);
	c[0] = dd_add_dd(d0, dd_mul_d(t, x));

	t = dd_add_dd(dd_mul_dpow2(d4, -4.0), dd_mul_d(d5, 5.0 * x));
	t = dd_add_dd(dd_mul_d(d3, 3.0), dd_mul_d(t, x));
	t = dd_add_dd(dd_mul_dpow2(d2, -2.0), dd_mul_d(t, x));
	c[1] = dd_add_dd(d1, dd_mul_d(t, x));

	t = dd_sub_dd(dd_mul_d(d4, 6.0), dd_mul_d(d5, 10.0 * x));
	t = dd_add_dd(dd_mul_d(d3, -3.0), dd_mul_d(t, x));
	c[2] = dd_add_dd(d2, dd_mul_d(t, x));

	t = dd_add_dd(dd_mul_dpow2(d4, -4.0), dd_mul_d(d5, 10.0 * x));
	c[3] = dd_add_dd(d3, dd_mul_d(t, x));

	c[4] = dd_add_dd(d4, dd_mul_d(d5, -5.0 * x));

	c[5] = d5;
}

/*-------------------------------------------------------------------------*/
static void
translate_gmp(curr_poly_t *c, mpz_t * gmp_z, 
		mpz_t * lin, int64 k)
{
	int64_2gmp(-k, c->gmp_help1);
	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[1], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[2], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
/* a0<-a0-a1*k+a2*k^2-a3*k^3+a4*k^4-a5*k^5 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[2], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)2);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)3);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)4);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)5);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
/* a1<-a1-2*a2*k+3*a3*k^2-4*a4*k^3+5*a5*k^4 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)3);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)6);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)10);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
/* a2<-a2-3*a3*k+6*a4*k^2-10*a5*k^3 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)4);
	mpz_add(gmp_z[3], gmp_z[3], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)10);
	mpz_add(gmp_z[3], gmp_z[3], c->gmp_help2);
/* a3<-a3-4*a4*k+10*a5*k^2 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, (mp_limb_t)5);
	mpz_add(gmp_z[4], gmp_z[4], c->gmp_help2);
/* a4<-a4-5*a5*k */

	mpz_mul(c->gmp_help2, lin[1], c->gmp_help1);
	mpz_add(lin[0], lin[0], c->gmp_help2);
/* lin0<-lin0-lin1*k */
}

/*-------------------------------------------------------------------------*/
#define TRANSLATE_SIZE 0
#define SKEWNESS 1

#define ROTATE0 0
#define ROTATE1 1
#define ROTATE2 2

typedef struct {
	ddpoly_t *rpoly;
	ddpoly_t *apoly;
	integrate_t *integ_aux;

	dickman_t *dickman_aux;
	double root_score_r;
	double root_score_a;

	double last_rotate_c0;
	double last_rotate_c1;
	double last_score;
} opt_data_t;

static double poly_xlate_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	double translated[MAX_POLY_DEGREE + 1];
	double t = floor(v[TRANSLATE_SIZE] + 0.5);
	double s = v[SKEWNESS];
	ddpoly_t *apoly = opt->apoly;

	translate_dbl(translated, apoly->coeff, t);

	return ifs(translated, s);
}

static double poly_rotate_callback(double *v, void *extra)
{
	opt_data_t *opt = (opt_data_t *)extra;
	ddpoly_t new_apoly = *(opt->apoly);
	dd_t r0 = opt->rpoly->coeff[0];
	dd_t r1 = opt->rpoly->coeff[1];
	double c0 = floor(v[ROTATE0] + 0.5);
	double c1 = floor(v[ROTATE1] + 0.5);
	double result = 0;

	if (c0 == opt->last_rotate_c0 &&
	    c1 == opt->last_rotate_c1) {
		return -(opt->last_score);
	}

	new_apoly.coeff[0] = dd_add_dd(new_apoly.coeff[0],
					dd_mul_d(r0, c0));
	new_apoly.coeff[1] = dd_add_dd(new_apoly.coeff[1],
					dd_add_dd(dd_mul_d(r0, c1),
						  dd_mul_d(r1, c0))
					);
	new_apoly.coeff[2] = dd_add_dd(new_apoly.coeff[2],
					dd_mul_d(r1, c1));

	analyze_poly_size(opt->integ_aux, opt->rpoly, &new_apoly, &result);

	opt->last_rotate_c0 = c0;
	opt->last_rotate_c1 = c1;
	opt->last_score = result;

	return -result;
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

	new_apoly.degree = 5;
	translate_dd(new_apoly.coeff, apoly->coeff, t);

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
	opt_data_t opt_data;
	int i;
	double best[MAX_VARS];
	double score;
	ddpoly_t rpoly, apoly;
	signed_mp_t tmp;

	opt_data.rpoly = &rpoly;
	opt_data.apoly = &apoly;
	opt_data.integ_aux = &curr_data->assess.integ_aux;
	opt_data.last_rotate_c0 = 1000000;
	opt_data.last_rotate_c1 = 0;

	rpoly.degree = 1;
	for (i = 0; i < 2; i++) {
		gmp2mp(c->gmp_lina[i], &tmp.num);
		tmp.sign = (mpz_sgn(c->gmp_lina[i]) < 0) ? NEGATIVE : POSITIVE;
		rpoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}
	apoly.degree = 5;
	for (i = 0; i < 6; i++) {
		gmp2mp(c->gmp_a[i], &tmp.num);
		tmp.sign = (mpz_sgn(c->gmp_a[i]) < 0) ? NEGATIVE : POSITIVE;
		apoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}

	best[ROTATE0] = 0;
	best[ROTATE1] = 0;

	score = minimize(best, 2, 1e-5, 40, poly_rotate_callback, &opt_data);

	mpz_set_d(c->gmp_help1, floor(best[ROTATE1] + 0.5));
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_a[1], c->gmp_a[1], c->gmp_help1);

	mpz_set_d(c->gmp_help1, floor(best[ROTATE0] + 0.5));
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[1], c->gmp_a[1], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_a[0], c->gmp_a[0], c->gmp_help1);

	*pol_norm = fabs(score);
#if 0
	printf("norm %.7e\n", *pol_norm);
	for (i = 0; i < 2; i++)
		gmp_printf("%+Zd\n", c->gmp_lina[i]);
	for (i = 0; i < 6; i++)
		gmp_printf("%+Zd\n", c->gmp_a[i]);
#endif
}

/*-------------------------------------------------------------------------*/
double
optimize_basic(dpoly_t *apoly, double *best_skewness,
		double *best_translation)
{
	uint32 i;
	opt_data_t opt_data;
	ddpoly_t new_apoly;
	double best[MAX_VARS];
	double score;
	double s0 = 1.0;
	double s1 = 1.0;

	new_apoly.degree = apoly->degree;
	for (i = 0; i <= apoly->degree; i++)
		new_apoly.coeff[i] = dd_set_d(apoly->coeff[i]);

	opt_data.rpoly = NULL;
	opt_data.apoly = &new_apoly;
	if (apoly->coeff[4] != 0)
		s0 = fabs(apoly->coeff[4] / apoly->coeff[5]);
	if (apoly->coeff[3] != 0)
		s1 = sqrt(fabs(apoly->coeff[3] / apoly->coeff[5]));

	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = MAX(s0, s1);

	score = minimize(best, 2, 1e-5, 40, poly_xlate_callback, &opt_data);

	*best_translation = floor(best[TRANSLATE_SIZE] + 0.5);
	*best_skewness  = best[SKEWNESS];
	return sqrt(score);
}

/*-------------------------------------------------------------------------*/
static void
optimize_final_core(curr_poly_t *c, assess_t *assess,
			double root_score, double *best_score_out,
			double *best_skewness_out)
{
	int i;
	signed_mp_t tmp;
	opt_data_t opt_data;
	double best[MAX_VARS];
	double score, curr_score;
	ddpoly_t rpoly, apoly;

	opt_data.rpoly = &rpoly;
	opt_data.apoly = &apoly;
	opt_data.integ_aux = &assess->integ_aux;
	opt_data.dickman_aux = &assess->dickman_aux;
	opt_data.root_score_r = 0.0;
	opt_data.root_score_a = root_score;
	rpoly.degree = 1;
	apoly.degree = 5;

	best[TRANSLATE_SIZE] = 0;
	best[SKEWNESS] = sqrt(fabs(mpz_get_d(c->gmp_b[3]) / 
				mpz_get_d(c->gmp_b[5])));
	best[SKEWNESS] = MAX(1.0, best[SKEWNESS]);
	curr_score = 0;

	do {
		score = curr_score;

		for (i = 0; i < 2; i++) {
			gmp2mp(c->gmp_linb[i], &tmp.num);
			tmp.sign = mpz_sgn(c->gmp_linb[i]) >= 0 ? 
						POSITIVE : NEGATIVE;
			rpoly.coeff[i] = dd_signed_mp2dd(&tmp);
		}
		for (i = 0; i < 6; i++) {
			gmp2mp(c->gmp_b[i], &tmp.num);
			tmp.sign = mpz_sgn(c->gmp_b[i]) >= 0 ? 
						POSITIVE : NEGATIVE;
			apoly.coeff[i] = dd_signed_mp2dd(&tmp);
		}

		curr_score = minimize(best, 2, 1e-5, 40, 
				poly_murphy_callback, &opt_data);

		translate_gmp(c, c->gmp_b, c->gmp_linb, 
				(int64)floor(best[TRANSLATE_SIZE] + 0.5));
		best[TRANSLATE_SIZE] = 0;

	} while (fabs(curr_score - score) > 0.05 * fabs(score));

	*best_score_out = fabs(curr_score);
	*best_skewness_out = best[SKEWNESS];
}

/*-------------------------------------------------------------------------*/
static void
get_bernstein_score(curr_poly_t *c, assess_t *assess,
		double root_score, double *eptr)
{
	int i;
	double size_score;
	ddpoly_t rpoly, apoly;
	signed_mp_t tmp;

	root_score = pow(exp(root_score), -2./6);

	rpoly.degree = 1;
	for (i = 0; i < 2; i++) {
		gmp2mp(c->gmp_linb[i], &tmp.num);
		tmp.sign = mpz_sgn(c->gmp_linb[i]) >= 0 ? POSITIVE : NEGATIVE;
		rpoly.coeff[i] = dd_signed_mp2dd(&tmp);
	}

	apoly.degree = 5;
	for (i = 0; i < 6; i++) {
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
	int i;
	double alpha, skewness, bscore, combined_score;
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	assess_t *assess = &s->assess;

	for (i = 0; i < 6; i++)
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
	for (i = 0; i < 2; i++)
		mpz_set(c->gmp_linb[i], c->gmp_lina[i]);

	if (stage2_root_score(5, c->gmp_b, data->murphy_p_bound, &alpha, 0))
		return;

	if (alpha > -4.5)
		return;

	get_bernstein_score(c, assess, alpha, &bscore);

#if 0
	printf("%.0lf %d %lf %le\n", (double)x, y, alpha, bscore); 
	fflush(stdout);
#endif

	if (bscore > data->min_e_bernstein) {

		optimize_final_core(c, assess, alpha, 
				&combined_score, &skewness);

		if (combined_score > data->min_e) {
			data->callback(data->callback_data, 5, c->gmp_b, 
					c->gmp_linb, skewness, bscore,
					alpha, combined_score);
		}
	}
}
