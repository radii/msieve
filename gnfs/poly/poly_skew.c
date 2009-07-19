/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 4/3/09
--------------------------------------------------------------------*/

#include "poly_skew.h"

#if MAX_POLY_DEGREE < 6
#error "Polynomial generation assumes degree <= 6 allowed"
#endif

typedef struct {
	double digits;
	double default_coeff;
	double stage1_norm;
	double stage2_norm;
	double final_norm;
	double bernstein_mult;
} poly_param_t;

typedef struct {
	poly_stage2_t *stage2_data;
	poly_config_t *config;
} stage1_callback_data_t;

	/* these parameters are shamelessly stolen from GGNFS;
	   the last item in the list is the minimum combined score
	   that polynomials must have, and we scale the GGNFS version
	   of the parameter because we use a different rating system */

static const poly_param_t prebuilt_params[] = {
	{ 98, 1,     2.60E+014, 1.40E+013, 3.45E-009, 2.20},
	{ 99, 1,     3.80E+014, 2.00E+013, 3.10E-009, 2.20},
	{100, 1,     5.56E+014, 2.86E+013, 2.80E-009, 2.20}, /* 3 ^ (.1 * x) */
	{101, 1,     8.13E+014, 4.08E+013, 2.50E-009, 2.46},
	{102, 1,     1.19E+015, 5.82E+013, 2.21E-009, 2.74},
	{103, 1,     1.74E+015, 8.31E+013, 1.94E-009, 3.06},
	{104, 1,     2.55E+015, 1.18E+014, 1.71E-009, 3.41},
	{105, 1,     3.74E+015, 1.69E+014, 1.50E-009, 3.81},
	{106, 1,     5.47E+015, 2.41E+014, 1.32E-009, 4.25},
	{107, 1,     8.01E+015, 3.44E+014, 1.16E-009, 4.75},
	{108, 1,     1.17E+016, 4.91E+014, 1.02E-009, 5.30},
	{109, 1,     1.72E+016, 7.00E+014, 8.92E-010, 5.91},
	{110, 1,     2.52E+016, 9.98E+014, 7.83E-010, 6.6}, /* 2.88 ^(.1*x) */
	{111, 1,     3.68E+016, 1.42E+015, 6.88E-010, 7.34},
	{112, 1,     5.39E+016, 2.03E+015, 6.04E-010, 8.15},
	{113, 1,     7.90E+016, 2.90E+015, 5.31E-010, 9.06},
	{114, 1,     1.16E+017, 4.13E+015, 4.66E-010, 10.1},
	{115, 1,     1.69E+017, 5.90E+015, 4.09E-010, 11.2},
	{116, 1,     2.48E+017, 8.41E+015, 3.60E-010, 12.4},
	{117, 1,     3.63E+017, 1.20E+016, 3.16E-010, 13.8},
	{118, 1,     5.31E+017, 1.71E+016, 2.77E-010, 15.4},
	{119, 1,     7.78E+017, 2.44E+016, 2.44E-010, 17.1},
	{120, 5000,  1.14E+018, 3.48E+016, 2.14E-010, 19.0}, /* 3 ^ (.1 * x) */
	{121, 8000,  1.67E+018, 4.97E+016, 1.88E-010, 21.2},
	{122, 10000, 2.44E+018, 7.09E+016, 1.65E-010, 23.7},
	{123, 12000, 3.58E+018, 1.01E+017, 1.45E-010, 26.4},
	{124, 14000, 5.24E+018, 1.44E+017, 1.27E-010, 29.5},
	{125, 16000, 7.67E+018, 2.06E+017, 1.12E-010, 32.9},
	{126, 18000, 1.12E+019, 2.94E+017, 9.82E-011, 36.7},
	{127, 20000, 1.64E+019, 4.19E+017, 8.62E-011, 41.0},
	{128, 23000, 2.41E+019, 5.98E+017, 7.57E-011, 45.8},
	{129, 26000, 3.52E+019, 8.52E+017, 6.65E-011, 51.1},
	{130, 29000, 5.16E+019, 1.22E+018, 5.84E-011, 57.0}, /* 3 ^ (.1*x)? */
	{131, 33000, 7.55E+019, 1.73E+018, 5.13E-011, 63.6},
	{132, 37000, 1.11E+020, 2.47E+018, 4.51E-011, 71.0},
	{133, 41000, 1.62E+020, 3.53E+018, 3.96E-011, 79.3},
	{134, 46000, 2.37E+020, 5.04E+018, 3.48E-011, 88.5},
	{135, 51000,3.47E+020, 7.18E+018, 3.05E-011, 98.7},
	{136, 56000,5.08E+020, 1.02E+019, 2.68E-011, 110.2},
	{137, 62000,7.44E+020, 1.46E+019, 2.36E-011, 123.0},
	{138, 68000,1.09E+021, 2.09E+019, 2.07E-011, 137.3},
	{139, 74000,1.60E+021, 2.98E+019, 1.82E-011, 153.2},
	{140, 81000,2.34E+021, 4.24E+019, 1.60E-011, 171.0},
	{141, 88000,3.42E+021, 6.05E+019, 1.40E-011, 190.8},
	{142, 95000,5.01E+021, 8.64E+019, 1.23E-011, 213.0},
	{143, 102000,7.33E+021, 1.23E+020, 1.08E-011, 237.7},
	{144, 110000,1.07E+022, 1.76E+020, 9.50E-012, 265.4},
	{145, 119000,1.57E+022, 2.51E+020, 8.34E-012, 296.2},
	{146, 129000,2.30E+022, 3.58E+020, 7.33E-012, 330.6},
	{147, 140000,3.37E+022, 5.10E+020, 6.43E-012, 369.0},
	{148, 152000,4.94E+022, 7.28E+020, 5.65E-012, 411.8},
	{149, 165000,7.23E+022, 1.04E+021, 4.96E-012, 459.6},
	{150, 179000,1.06E+023, 1.48E+021, 4.36E-012, 513.0},
	{151, 194000,1.55E+023, 2.11E+021, 3.83E-012, 572.6},
	{152, 210000,2.27E+023, 3.01E+021, 3.36E-012, 639.1},
	{153, 227000,3.32E+023, 4.30E+021, 2.95E-012, 713.3},
	{154, 245000,4.86E+023, 6.13E+021, 2.59E-012, 796.1},
	{155, 264000,7.12E+023, 8.75E+021, 2.28E-012, 888.5},
};

/*--------------------------------------------------------------------*/
static void get_poly_params(double digits, poly_param_t *params) {

	const poly_param_t *low, *high;
	uint32 i, max_entry;
	double j, k, dist;
	double max_digits;

	/* if the input is too small (large), just use 
	   the first (last) table entry */

	if (digits < prebuilt_params[0].digits) {
		*params = prebuilt_params[0];
		return;
	}

	max_entry = sizeof(prebuilt_params) / sizeof(poly_param_t);
	max_digits = prebuilt_params[max_entry - 1].digits;
	if (digits >= max_digits) {
		if (digits > max_digits + 5) {
			printf("error: no parameters for "
				"%.0lf digit inputs\n", digits + 0.5);
			exit(-1);
		}
		*params = prebuilt_params[max_entry - 1];
		return;
	}

	/* Otherwise the parameters to use are a weighted average 
	   of the two table entries the input falls between */

	for (i = 0; i < max_entry - 1; i++) {
		if (digits < prebuilt_params[i+1].digits)
			break;
	}

	low = &prebuilt_params[i];
	high = &prebuilt_params[i+1];
	dist = high->digits - low->digits;
	j = digits - low->digits;
	k = high->digits - digits;

	params->digits = digits;
	params->default_coeff = (low->default_coeff * k +
			       high->default_coeff * j) / dist;
	params->stage1_norm = (low->stage1_norm * k +
			       high->stage1_norm * j) / dist;
	params->stage2_norm = (low->stage2_norm * k +
			       high->stage2_norm * j) / dist;
	params->final_norm = (low->final_norm * k +
			       high->final_norm * j) / dist;
	params->bernstein_mult = (low->bernstein_mult * k +
			       high->bernstein_mult * j) / dist;
}

/*------------------------------------------------------------------*/
static void stage1_callback(mpz_t a5, mpz_t p, mpz_t m, 
				double a3_bound, void *extra) {
	
	stage1_callback_data_t *d = (stage1_callback_data_t *)extra;

	poly_stage2_run(d->stage2_data, a5, p, m, a3_bound);
}

/*------------------------------------------------------------------*/
static void
stage2_callback(void *extra, int deg, 
		mpz_t * coeff1, mpz_t * coeff2,
		double skewness, double size_score,
		double root_score, double combined_score)
{
	int i;
	poly_select_t poly;
	mp_poly_t *rpoly;
	mp_poly_t *apoly;
	poly_config_t *config = (poly_config_t *)extra;

	memset(&poly, 0, sizeof(poly_select_t));
	rpoly = &poly.rpoly;
	apoly = &poly.apoly;

	rpoly->degree = 1;
	for (i = 0; i <= 1; i++) {
		gmp2mp(coeff2[i], &rpoly->coeff[i].num);
		rpoly->coeff[i].sign = POSITIVE;
		if (mpz_sgn(coeff2[i]) < 0)
			rpoly->coeff[i].sign = NEGATIVE;
	}

	apoly->degree = deg;
	for (i = 0; i <= deg; i++) {
		gmp2mp(coeff1[i], &apoly->coeff[i].num);
		apoly->coeff[i].sign = POSITIVE;
		if (mpz_sgn(coeff1[i]) < 0)
			apoly->coeff[i].sign = NEGATIVE;
	}
	poly.root_score = root_score;
	poly.size_score = size_score;
	poly.combined_score = combined_score;
	poly.skewness = skewness;

	printf("save %le %lf %lf %le\n", size_score,
			root_score, skewness, combined_score);
	fprintf(config->all_poly_file, 
		"# norm %le alpha %lf e %.3le\nskew: %.2lf\n", 
		size_score, root_score, combined_score, skewness);
	for (i = 0; i <= deg; i++) {
		fprintf(config->all_poly_file, "c%u: %s", i,
				mpz_sgn(coeff1[i]) >= 0 ? " " : "");
		mpz_out_str(config->all_poly_file, 10, coeff1[i]);
		fprintf(config->all_poly_file, "\n");
	}
	for (i = 0; i <= 1; i++) {
		fprintf(config->all_poly_file, "Y%u: %s", i,
				mpz_sgn(coeff2[i]) >= 0 ? " " : "");
		mpz_out_str(config->all_poly_file, 10, coeff2[i]);
		fprintf(config->all_poly_file, "\n");
	}
	fflush(config->all_poly_file);

	save_poly(config, &poly);
}

/*------------------------------------------------------------------*/
static void find_poly_core(msieve_obj *obj, mp_t *n,
			poly_config_t *config,
			uint32 degree, uint32 deadline) {

	double dbl_N, digits;
	poly_stage1_t stage1_data;
	poly_stage2_t stage2_data;
	poly_param_t params;
	stage1_callback_data_t stage1_callback_data;

	degree = 5;    /* required */

	/* initialize structures */

	stage1_callback_data.stage2_data = &stage2_data;
	stage1_callback_data.config = config;

	poly_stage1_init(&stage1_data, stage1_callback, 
			&stage1_callback_data);

	poly_stage2_init(&stage2_data, stage2_callback, config);

	/* get poly selection parameters */

	mp2gmp(n, stage1_data.gmp_N);
	dbl_N = mpz_get_d(stage1_data.gmp_N);
	digits = log(dbl_N) / log(10.0);
	get_poly_params(digits, &params);

	/* fill stage 1 data */

	stage1_data.norm_max = params.stage1_norm;
	stage1_data.deadline = deadline;

	if (obj->nfs_lower)
		uint64_2gmp(obj->nfs_lower, stage1_data.gmp_a5_begin);
	else
		mpz_set_d(stage1_data.gmp_a5_begin, params.default_coeff);

	if (obj->nfs_upper)
		uint64_2gmp(obj->nfs_upper, stage1_data.gmp_a5_end);
	else
		mpz_set_d(stage1_data.gmp_a5_end,
			pow(dbl_N, 1.0 / (double)(degree * 
					(degree - 1))) / 30 );

	if (obj->nfs_lower && obj->nfs_upper)
		stage1_data.deadline = 0;

	logprintf(obj, "searching leading coefficients from %.0lf to %.0lf\n",
			mpz_get_d(stage1_data.gmp_a5_begin),
			mpz_get_d(stage1_data.gmp_a5_end));

	/* fill stage 2 data */

	mp2gmp(n, stage2_data.gmp_N);
	stage2_data.max_norm = params.stage2_norm;
	stage2_data.min_e = params.final_norm;
	stage2_data.min_e_bernstein = params.final_norm / 
					params.bernstein_mult;

	poly_stage1_run(obj, &stage1_data);

	poly_stage1_free(&stage1_data);
	poly_stage2_free(&stage2_data);
}

/*------------------------------------------------------------------*/
void find_poly_skew(msieve_obj *obj, mp_t *n,
			poly_config_t *config,
			uint32 deadline) {

	/* search for NFS polynomials */

	uint32 bits = mp_bits(n);

	if (bits < 331) {		/* <= 100 digits */
		find_poly_core(obj, n, config, 4, deadline);
	}
	else if (bits < 347) {		/* 100-105 digits */
//		find_poly_core(obj, n, config, 4, deadline);
		find_poly_core(obj, n, config, 5, deadline);
	}
	else if (bits < 592) {		/* 105-180 digits */
		find_poly_core(obj, n, config, 5, deadline);
	}
	else if (bits < 659) {		/* 180-200 digits */
		find_poly_core(obj, n, config, 5, deadline);
//		find_poly_core(obj, n, config, 6, deadline);
	}
	else {				/* 200+ digits */
		find_poly_core(obj, n, config, 6, deadline);
	}
}
