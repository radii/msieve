#include "stage1.h"

/*------------------------------------------------------------------------*/
void
p_batch_init(p_batch_t *pb) {

	pb->num_entries = 0;
	pb->num_small_entries = 0;
	pb->num_entries_alloc = 1000;
	pb->entries = (p_entry_t *)xmalloc(pb->num_entries_alloc *
					sizeof(p_entry_t));
	pb->num_index = 0;
	pb->num_index_alloc = 1000;
	pb->index = (index_t *)xmalloc(pb->num_index_alloc *
					sizeof(index_t));
	pb->num_roots = 0;
	pb->num_roots_alloc = 5000;
	pb->roots = (uint64 *)xmalloc(pb->num_roots_alloc *
					sizeof(uint64));
}

/*------------------------------------------------------------------------*/
void
p_batch_free(p_batch_t *pb) 
{
	free(pb->entries);
	free(pb->index);
	free(pb->roots);
}

/*------------------------------------------------------------------------*/
void
p_batch_reset(p_batch_t *pb) 
{
	pb->num_entries = 0;
	pb->num_small_entries = 0;
	pb->num_index = 0;
	pb->num_roots = 0;
}

/*------------------------------------------------------------------------*/
void
p_batch_add(p_batch_t *pb, uint64 p, uint32 poly_idx,
		uint32 num_roots, mpz_t *roots) 
{
	uint32 i;
	p_entry_t *entry;
	index_t *index;
	uint64 *r;

	if (pb->num_entries == 0 ||
	    pb->entries[pb->num_entries - 1].p != p) {

		if (pb->num_entries == pb->num_entries_alloc) {
			pb->num_entries_alloc *= 2;
			pb->entries = (p_entry_t *)xrealloc(pb->entries,
						pb->num_entries_alloc *
						sizeof(p_entry_t));
		}
		entry = pb->entries + pb->num_entries++;
		entry->p = p;
		entry->num_poly = 1;
		entry->num_roots = num_roots;
		entry->index_start_offset = pb->num_index;
		if (p < 65536)
			pb->num_small_entries++;
	}
	else {
		entry = pb->entries + pb->num_entries - 1;
		entry->num_poly++;

		if (num_roots != entry->num_roots) {
			printf("error: mismatch in number of roots\n");
			exit(-1);
		}
	}

	if (pb->num_index == pb->num_index_alloc) {
		pb->num_index_alloc *= 2;
		pb->index = (index_t *)xrealloc(pb->index,
						pb->num_index_alloc *
						sizeof(index_t));
	}
	index = pb->index + pb->num_index++;
	index->which_poly = poly_idx;
	index->start_offset = pb->num_roots;

	if (pb->num_roots + num_roots >= pb->num_roots_alloc) {
		pb->num_roots_alloc = 2 * (pb->num_roots + num_roots);
		pb->roots = (uint64 *)xrealloc(pb->roots,
						pb->num_roots_alloc *
						sizeof(uint64));
	}
	r = pb->roots + pb->num_roots;
	for (i = 0; i < num_roots; i++)
		r[i] = gmp2uint64(roots[i]);
	pb->num_roots += num_roots;
}

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
poly_batch_init(poly_batch_t *poly, poly_stage1_t *data)
{
	uint32 i;

	for (i = 0; i < POLY_BATCH_SIZE; i++) {
		mpz_init(poly->batch[i].high_coeff);
		mpz_init(poly->batch[i].trans_N);
		mpz_init(poly->batch[i].trans_m0);
	}

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
poly_batch_free(poly_batch_t *poly)
{
	uint32 i;

	for (i = 0; i < POLY_BATCH_SIZE; i++) {
		mpz_clear(poly->batch[i].high_coeff);
		mpz_clear(poly->batch[i].trans_N);
		mpz_clear(poly->batch[i].trans_m0);
	}

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
search_coeffs_core(msieve_obj *obj, poly_batch_t *poly, uint32 digits)
{
	uint32 i, j;
	uint32 deadline;
	uint32 degree = poly->degree;
	uint32 mult = 0;

	switch (degree) {
	case 4: mult = 4 * 4 * 4 * 4; break;
	case 5: mult = 5 * 5 * 5 * 5 * 5; break;
	case 6: mult = 6 * 6 * 6 * 6 * 6 * 6; break;
	}

	for (i = 0; i < poly->num_poly; i++) {
		curr_poly_t *c = poly->batch + i;
		mpz_mul_ui(c->trans_N, poly->N, (mp_limb_t)mult);
		for (j = 0; j < degree - 1; j++)
			mpz_mul(c->trans_N, c->trans_N, c->high_coeff);

		mpz_root(c->trans_m0, c->trans_N, (mp_limb_t)degree);
		mpz_tdiv_q(poly->m0, poly->N, c->high_coeff);
		mpz_root(poly->m0, poly->m0, (mp_limb_t)degree);

		c->sieve_size = c->coeff_max / mpz_get_d(poly->m0) * 
				c->p_size_max * c->p_size_max / degree;
	}

	deadline = 800;
	if (digits <= 100)
		deadline = 5;
	else if (digits <= 105)
		deadline = 20;
	else if (digits <= 110)
		deadline = 30;
	else if (digits <= 120)
		deadline = 50;
	else if (digits <= 130)
		deadline = 100;
	else if (digits <= 140)
		deadline = 200;
	else if (digits <= 150)
		deadline = 400;
	printf("deadline: %u seconds per coefficient\n", deadline);

	sieve_lattice(obj, poly, 2000, 2001, 100000, deadline * i);
}

/*------------------------------------------------------------------------*/
static void
search_coeffs(msieve_obj *obj, poly_batch_t *poly, 
		bounds_t *bounds, uint32 deadline)
{
	mpz_t curr_high_coeff;
	double dn = mpz_get_d(poly->N);
	uint32 digits = mpz_sizeinbase(poly->N, 10);
	uint32 batch_size = POLY_BATCH_SIZE;
	double start_time = get_cpu_time();

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

	if (mpz_cmp_ui(curr_high_coeff, (mp_limb_t)2000) < 0)
		batch_size = MIN(batch_size, 10);

	poly->num_poly = 0;
	while (1) {
		curr_poly_t *c = poly->batch + poly->num_poly;

		if (mpz_cmp(curr_high_coeff, bounds->gmp_high_coeff_end) > 0) {
			if (poly->num_poly)
				search_coeffs_core(obj, poly, digits);
			break;
		}

		stage1_bounds_update(bounds, dn, 
					mpz_get_d(curr_high_coeff),
					poly->degree);
		mpz_set(c->high_coeff, curr_high_coeff);
		c->p_size_max = bounds->p_size_max;
		c->coeff_max = bounds->coeff_max;

		if (++poly->num_poly == batch_size) {
			search_coeffs_core(obj, poly, digits);

			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;

			if (deadline) {
				double curr_time = get_cpu_time();
				double elapsed = curr_time - start_time;

				if (elapsed > deadline)
					break;
			}

			poly->num_poly = 0;
			if (mpz_cmp_ui(curr_high_coeff, (mp_limb_t)2000) >= 0)
				batch_size = POLY_BATCH_SIZE;
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
	poly_batch_t poly;

	stage1_bounds_init(&bounds, data);
	poly_batch_init(&poly, data);

	search_coeffs(obj, &poly, &bounds, data->deadline);

	stage1_bounds_free(&bounds);
	poly_batch_free(&poly);
	return 1;
}
