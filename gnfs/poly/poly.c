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

#include "poly.h"

/* used to place a deadline on how long polynomial 
   selection will run. Note that the time budget is
   independent of CPU speed; faster CPUs will simply
   search more of the polynomial space */

typedef struct {
	uint32 bits;
	uint32 seconds;
} poly_deadline_t;

static const poly_deadline_t time_limits[] = {
	{MIN_NFS_BITS, 4 * 60},
	{304, 8 * 60},
	{320, 15 * 60},
	{348, 30 * 60},
	{365, 1 * 3600},
	{383, 2 * 3600},
	{399, 4 * 3600},
	{416, 8 * 3600},
	{433, 16 * 3600},
	{449, 32 * 3600},
	{466, 64 * 3600},
	{482, 100 * 3600},
	{498, 200 * 3600},
	{514, 300 * 3600},
};

#define NUM_TIME_LIMITS sizeof(time_limits)/sizeof(time_limits[0])

/*------------------------------------------------------------------*/
int32 read_poly(msieve_obj *obj, mp_t *n,
	       mp_poly_t *rat_poly,
	       mp_poly_t *alg_poly,
	       double *skewness) {
	
	int32 i;
	FILE *fp;
	char buf[LINE_BUF_SIZE];
	mp_t read_n;
	signed_mp_t val, rpow, tmp;

	fp = fopen(obj->nfs_fbfile_name, "r");
	if (fp == NULL)
		return -1;
	
	buf[0] = 0;
	fgets(buf, (int)sizeof(buf), fp);
	if (buf[0] != 'N') {
		fclose(fp);
		logprintf(obj, "warning: factor base file uninitialized\n");
		return -1;
	}

	/* check that the polynomial is for the 
	   right number */

	mp_str2mp(buf + 2, &read_n, 10);
	if (mp_cmp(&read_n, n)) {
		fclose(fp);
		logprintf(obj, "warning: NFS input not found in "
				"factor base file\n");
		return -1;
	}

	/* read in skewness if present */

	fgets(buf, (int)sizeof(buf), fp);
	if (buf[0] == 'S') {
		if (skewness != NULL)
			*skewness = atof(buf + 5);
		fgets(buf, (int)sizeof(buf), fp);
	}
	else if (skewness != NULL) {
		*skewness = 1.0;
	}

	/* read one coefficient per line; 'R<number>' is
	   for rational coefficients, 'A<number>' for algebraic */

	while ((buf[0] == 'R' || buf[0] == 'A') && isdigit(buf[1])) {
		signed_mp_t *read_coeff;
		char *tmp;

		i = buf[1] - '0';
		if (i > MAX_POLY_DEGREE) {
			fclose(fp);
			logprintf(obj, "warning: polynomial degree exceeds "
					"%d\n", MAX_POLY_DEGREE);
			exit(-1);
		}

		if (buf[0] == 'R')
			read_coeff = rat_poly->coeff + i;
		else
			read_coeff = alg_poly->coeff + i;

		tmp = buf + 2;
		while (isspace(*tmp))
			tmp++;

		if (*tmp == '-') {
			read_coeff->sign = NEGATIVE;
			tmp++;
		}
		else {
			read_coeff->sign = POSITIVE;
		}
		mp_str2mp(tmp, &read_coeff->num, 10);
		if (fgets(buf, (int)sizeof(buf), fp) == NULL)
			break;
	}

	for (i = MAX_POLY_DEGREE; i >= 0; i--) {
		if (!mp_is_zero(&rat_poly->coeff[i].num))
			break;
	}
	if (i > 0)
		rat_poly->degree = i;

	for (i = MAX_POLY_DEGREE; i >= 0; i--) {
		if (!mp_is_zero(&alg_poly->coeff[i].num))
			break;
	}
	if (i > 0)
		alg_poly->degree = i;

	fclose(fp);

	if (rat_poly->degree == 0 || alg_poly->degree == 0) {
		logprintf(obj, "error: polynomial is missing or corrupt\n");
		exit(-1);
	}
	if (rat_poly->degree != 1) {
		logprintf(obj, "error: no support for nonlinear "
				"rational polynomials\n");
		exit(-1);
	}
	
	/* plug the rational polynomial coefficients into the 
	   algebraic polynomial */

	i = alg_poly->degree;
	signed_mp_copy(alg_poly->coeff + i, &val);
	signed_mp_copy(rat_poly->coeff + 1, &rpow);

	for (i--; i >= 0; i--) {
		signed_mp_mul(&val, rat_poly->coeff + 0, &tmp);
		tmp.sign = (tmp.sign == POSITIVE) ? NEGATIVE : POSITIVE;
		signed_mp_copy(&tmp, &val);

		signed_mp_mul(alg_poly->coeff + i, &rpow, &tmp);
		signed_mp_add(&val, &tmp, &val);

		signed_mp_mul(rat_poly->coeff + 1, &rpow, &tmp);
		signed_mp_copy(&tmp, &rpow);
	}

	/* verify that |result| >= N, and that result % N == 0. 
	   The only place where we do any mod-N arithmetic is the 
	   NFS square root, which will not work if N has additional 
	   factors that are not reflected in the polynomials */

	if ((i = mp_cmp(&val.num, n)) < 0) {
		logprintf(obj, "error: NFS input does not match polynomials\n");
		logprintf(obj, "check that input doesn't have small factors\n");
		exit(-1);
	}
	else if (i > 0) {
		mp_mod(&val.num, n, &read_n);
		if (!mp_is_zero(&read_n)) {
			logprintf(obj, "error: NFS input does not "
					"match polynomials\n");
			exit(-1);
		}
	}

	return 0;
}

/*------------------------------------------------------------------*/
void write_poly(msieve_obj *obj, mp_t *n,
	       mp_poly_t *rat_poly,
	       mp_poly_t *alg_poly,
	       double skewness) {
	
	/* log a generated polynomial to the factor base file */

	uint32 i;
	FILE *fp;

	fp = fopen(obj->nfs_fbfile_name, "w");
	if (fp == NULL) {
		printf("error; cannot open factor base file '%s'\n",
					obj->nfs_fbfile_name);
		exit(-1);
	}

	fprintf(fp, "N %s\n", mp_sprintf(n, 10, obj->mp_sprintf_buf));
	if (skewness > 0)
		fprintf(fp, "SKEW %.2lf\n", skewness);

	for (i = 0; i <= rat_poly->degree; i++) {
		fprintf(fp, "R%u %c%s\n", i,
			(rat_poly->coeff[i].sign == POSITIVE? ' ' : '-'),
			mp_sprintf(&rat_poly->coeff[i].num, 10,
						obj->mp_sprintf_buf));
	}
	for (i = 0; i <= alg_poly->degree; i++) {
		fprintf(fp, "A%u %c%s\n", i,
			(alg_poly->coeff[i].sign == POSITIVE? ' ' : '-'),
			mp_sprintf(&alg_poly->coeff[i].num, 10,
						obj->mp_sprintf_buf));
	}
	fclose(fp);
}

/*------------------------------------------------------------------*/
int32 find_poly(msieve_obj *obj, mp_t *n) {

	/* external entry point for NFS polynomial generation */

	uint32 i, j;
	poly_config_t config;
	uint32 deadline;

	logprintf(obj, "commencing number field sieve polynomial selection\n");

	/* do sanity checking */

	if ((obj->nfs_lower == 0 && obj->nfs_upper != 0) ||
	    (obj->nfs_lower != 0 && obj->nfs_upper == 0) ) {
		printf("lower/upper bounds must both be specified\n");
		return -3;
	}

	poly_config_init(&config);

	/* figure out how long poly selection should take */

	i = mp_bits(n);
	for (j = 0; j < NUM_TIME_LIMITS; j++) {
		if (i < time_limits[j].bits)
			break;
	}
	if (j == NUM_TIME_LIMITS) {
		deadline = time_limits[j-1].seconds;
	}
	else {
		const poly_deadline_t *low = &time_limits[j-1];
		const poly_deadline_t *high = &time_limits[j];
		uint32 dist = high->bits - low->bits;
		deadline = (uint32)(
			 ((double)low->seconds * (high->bits - i) +
			  (double)high->seconds * (i - low->bits)) / dist);
	}

	/* run the core polynomial finder */

	obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	find_poly_skew(obj, n, &config, deadline);
	obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	/* save the best polynomial */

	logprintf(obj, "polynomial selection complete\n");
	if (config.heap[0]->rpoly.degree > 0) {
		write_poly(obj, n, &config.heap[0]->rpoly,
				&config.heap[0]->apoly,
				config.heap[0]->skewness);
	}
	poly_config_free(&config);

	return 0;
}
