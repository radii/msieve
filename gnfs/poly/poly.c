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

	fp = fopen(obj->nfs_fbfile_name, "r");
	if (fp == NULL)
		return -1;
	
	/* check that the polynomial is for the 
	   right number */

	fgets(buf, (int)sizeof(buf), fp);
	if (buf[0] != 'N') {
		fclose(fp);
		return -1;
	}
	mp_str2mp(buf + 2, &read_n, 10);
	if (mp_cmp(&read_n, n)) {
		fclose(fp);
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

	while (!feof(fp) && (buf[0] == 'R' || buf[0] == 'A')) {
		signed_mp_t *read_coeff;
		char *tmp;

		if (buf[0] == 'R')
			read_coeff = rat_poly->coeff + (buf[1] - '0');
		else
			read_coeff = alg_poly->coeff + (buf[1] - '0');

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
		fgets(buf, (int)sizeof(buf), fp);
	}

	/* do some sanity checking */

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
		logprintf(obj, "warning: polynomial is corrupt\n");
		return -1;
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
	logprintf(obj, "time limit set to %.2f hours\n", deadline / 3600.0);

	/* run the core polynomial finder */

	obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
#if defined(HAVE_GMP)
	find_poly_skew(obj, n, &config, deadline);
#else
	find_poly_noskew(obj, n, &config, deadline);
#endif
	obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	/* save the best polynomial */

	logprintf(obj, "polynomial selection complete\n");
	write_poly(obj, n, &config.heap[0]->rpoly,
			&config.heap[0]->apoly,
			config.heap[0]->skewness);
	poly_config_free(&config);

	return 0;
}
