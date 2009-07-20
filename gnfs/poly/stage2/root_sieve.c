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

#define ROOT_HEAP_SIZE 1000
#define LATTICE_HEAP_SIZE 20
#define MAX_SIEVE_PRIME 100
#define MAX_SIEVE_PRIME_POWER 1000
#define LOG_SCALE_FACTOR 1000
#define UNROLL 4
#define DEFAULT_BLOCK_SIZE  8192
#define ROOT_SCORE_COARSE_MIN (-4.0)

static const double good_alpha[3][3] = {
	{-4.4, -4.8, -5.2},
	{-4.8, -5.6, -6.3},
	{-5.6, -6.2, -7.0},
};

/*-------------------------------------------------------------------------*/
static double
init_sieve_core(curr_poly_t *c, sieve_power_t *sp, uint32 deg, uint32 p)
{
	uint32 i, j;
	uint32 fi, mi, md, mp, step;
	uint32 co[MAX_POLY_DEGREE + 1];
	double contrib = 0;
	uint32 num_roots = 0;
	uint32 pk = sp->power;

	for (i = 0; i <= deg; i++)
		co[i] = mpz_fdiv_ui(c->gmp_a[i], (mp_limb_t)pk);
	md = mpz_fdiv_ui(c->gmp_d, (mp_limb_t)pk);
	mp = mpz_fdiv_ui(c->gmp_p, (mp_limb_t)pk);

	for (i = 0; i < pk; i++) {
		mi = mp_modsub_1(md, mp_modmul_1(i, mp, pk), pk);

		fi = co[deg];
		for (j = deg; j; j--)
			fi = (fi * i + co[j - 1]) % pk;

		if (mi != 0) {
			step = pk;
			while (mi % p == 0 && fi % p == 0) {
				mi /= p;
				fi /= p;
				step /= p;
			}
			if (mi % p) {
				sieve_root_t *r = sp->roots + num_roots++;
				r->resclass = i;
				r->step = step;
				r->start = mp_modmul_1(fi, 
						mp_modinv_1(mi, step), step);
			}	/* otherwise no solutions */
		}
		else if (fi == 0) {
			/* p^k | fi+j*(p*i-m) for all j */
			contrib += sp->root_contrib;
		}
	}

	sp->num_roots = num_roots;
	return contrib;
}

/*-------------------------------------------------------------------------*/
static void
init_sieve(curr_poly_t *c, root_sieve_t *rs, 
		uint32 deg, double sieve_bias)
{
	uint32 i, j;

	for (i = 0; i < rs->num_primes; i++) {
		sieve_prime_t *sp = rs->primes + i;

		for (j = 0; j < sp->num_powers; j++) {
			sieve_bias += init_sieve_core(c, sp->powers + j,
							deg, sp->prime);
		}
	}

	rs->sieve_bias = sieve_bias;
}

/*-------------------------------------------------------------------------*/
static void
compute_line_size(double max_norm, double *dbl_a, uint32 deg,
		  double dbl_p, double dbl_d, 
		  int64 last_xmin_in, int64 last_xmax_in,
		  int64 *xmin, int64 *xmax)
{
	uint32 i;
	double v0;
	dpoly_t apoly;
	double new_xlate, new_skewness;
	double x0, x1, offset;
	double last_xmin = (double)last_xmin_in;
	double last_xmax = (double)last_xmax_in;

	apoly.degree = deg;
	for (i = 0; i <= deg; i++)
		apoly.coeff[i] = dbl_a[i];

	apoly.coeff[1] = dbl_a[1] + dbl_p * last_xmin;
	apoly.coeff[0] = dbl_a[0] - dbl_d * last_xmin;
	v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
	offset = 10000;
	if (v0 > max_norm) {
		x0 = last_xmin;
		x1 = last_xmin + offset;
		while (1) {
			double new_x1 = last_xmin + offset;
			apoly.coeff[1] = dbl_a[1] + dbl_p * new_x1;
			apoly.coeff[0] = dbl_a[0] - dbl_d * new_x1;
			v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
			if (v0 <= max_norm || new_x1 >= last_xmax) {
				x1 = new_x1;
				break;
			}
			x0 = new_x1;
			offset *= 2;
		}
	}
	else {
		x0 = last_xmin - offset;
		x1 = last_xmin;
		while (1) {
			double new_x0 = last_xmin - offset;
			apoly.coeff[1] = dbl_a[1] + dbl_p * new_x0;
			apoly.coeff[0] = dbl_a[0] - dbl_d * new_x0;
			v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
			if (v0 > max_norm) {
				x0 = new_x0;
				break;
			}
			x1 = new_x0;
			offset *= 2;
		}
	}

	while (x1 - x0 > 500) {
		double xx = (x0 + x1) / 2;
		apoly.coeff[1] = dbl_a[1] + dbl_p * xx;
		apoly.coeff[0] = dbl_a[0] - dbl_d * xx;
		v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
		if (v0 > max_norm)
			x0 = xx;
		else
			x1 = xx;
	}
	*xmin = (int64)x0;

	apoly.coeff[1] = dbl_a[1] + dbl_p * last_xmax;
	apoly.coeff[0] = dbl_a[0] - dbl_d * last_xmax;
	v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
	offset = 10000;
	if (v0 > max_norm) {
		x0 = last_xmax - offset;
		x1 = last_xmax;
		while (1) {
			double new_x0 = last_xmax - offset;
			apoly.coeff[1] = dbl_a[1] + dbl_p * new_x0;
			apoly.coeff[0] = dbl_a[0] - dbl_d * new_x0;
			v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
			if (v0 <= max_norm || new_x0 <= last_xmin) {
				x0 = new_x0;
				break;
			}
			x1 = new_x0;
			offset *= 2;
		}
	}
	else {
		x0 = last_xmax;
		x1 = last_xmax + offset;
		while (1) {
			double new_x1 = last_xmax + offset;
			apoly.coeff[1] = dbl_a[1] + dbl_p * new_x1;
			apoly.coeff[0] = dbl_a[0] - dbl_d * new_x1;
			v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
			if (v0 > max_norm) {
				x1 = new_x1;
				break;
			}
			x0 = new_x1;
			offset *= 2;
		}
	}

	while (x1 - x0 > 500) {
		double xx = (x0 + x1) / 2;
		apoly.coeff[1] = dbl_a[1] + dbl_p * xx;
		apoly.coeff[0] = dbl_a[0] - dbl_d * xx;
		v0 = optimize_basic(&apoly, &new_skewness, &new_xlate);
		if (v0 > max_norm)
			x1 = xx;
		else
			x0 = xx;
	}
	*xmax = (double)x1;
}

/*-------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify_rotation(root_heap_t *heap, uint32 index, uint32 size) {

	uint32 c;
	rotation_t tmp;
	rotation_t *h = heap->entries;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].score > h[c+1].score)
			c++;

		if (h[index].score > h[c].score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].score > h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap_rotation(root_heap_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify_rotation(h, (uint32)i, size);
}

static void
save_rotation(root_heap_t *heap, int64 x, int32 y, uint32 score)
{
	rotation_t *h = heap->entries;

	if (heap->extra != NULL) {

		int32 abs_y = abs(y);
		int64 abs_x = x;

		if (x < 0)
			abs_x = -abs_x;

		if (abs_y > heap->cutoffs[1].y ||
		    abs_x > heap->cutoffs[1].x) {
			if (score > heap->default_cutoff) {
				optimize_final(x, y, 
					(poly_stage2_t *)heap->extra);
				return;
			}
		}
		else if (abs_y > heap->cutoffs[0].y ||
		         abs_x > heap->cutoffs[0].x) {
			score += heap->default_cutoff -
				  heap->cutoffs[1].score;
		}
		else {
			score += heap->default_cutoff -
				  heap->cutoffs[0].score;
		}

		if (score > heap->default_cutoff) {
			optimize_final(x, y, (poly_stage2_t *)heap->extra);
			return;
		}
	}

	if (heap->num_entries <= heap->max_entries - 1) {
		rotation_t *r = h + heap->num_entries++;
		r->x = x;
		r->y = y;
		r->score = score;
		if (heap->num_entries == heap->max_entries)
			make_heap_rotation(heap, heap->max_entries);
	}
	else if (h[0].score < score) {
		h[0].x = x;
		h[0].y = y;
		h[0].score = score;
		heapify_rotation(heap, 0, heap->max_entries);
	}
}

/*-------------------------------------------------------------------------*/
static void
sieve_one_block(uint16 *sieve_block, uint32 sieve_block_size,
		sieve_prime_t *primes, uint32 num_primes)
{
	uint32 i, j;
	uint64 *block_array = (uint64 *)sieve_block;
	uint32 block_words = sieve_block_size / UNROLL;

	memset(block_array, 0, block_words * sizeof(uint64));

	for (i = 0; i < num_primes; i++) {
		sieve_prime_t *sp = primes + i;
		uint64 *contrib_array = (uint64 *)sp->contrib_array;
		uint32 contrib_words = sp->contrib_array_size;
		uint32 contrib_offset = sp->contrib_array_offset;
		uint32 block_offset = 0;

		while (block_offset < block_words) {
			uint32 curr_words = MIN(contrib_words - contrib_offset,
						block_words - block_offset);
			uint64 *b = block_array + block_offset;
			uint64 *c = contrib_array + contrib_offset;

#if defined(GCC_ASM32X) && defined(HAS_MMX)
			j = 0;
			ASM_G volatile(
			    "cmpl $0, %3              \n\t"
			    "je 1f                   \n\t"
			    ALIGN_LOOP
			    "0:                       \n\t"
			    "movq (%1,%0,8), %%mm0    \n\t"
			    "movq 8(%1,%0,8), %%mm1   \n\t"
			    "movq 16(%1,%0,8), %%mm2  \n\t"
			    "movq 24(%1,%0,8), %%mm3  \n\t"
			    "paddw (%2,%0,8), %%mm0   \n\t"
			    "paddw 8(%2,%0,8), %%mm1  \n\t"
			    "paddw 16(%2,%0,8), %%mm2 \n\t"
			    "paddw 24(%2,%0,8), %%mm3 \n\t"
			    "movq %%mm0, (%1,%0,8)    \n\t"
			    "movq %%mm1, 8(%1,%0,8)   \n\t"
			    "movq %%mm2, 16(%1,%0,8)  \n\t"
			    "movq %%mm3, 24(%1,%0,8)  \n\t"
			    "addl $4, %0              \n\t"
			    "cmpl %3, %0              \n\t"
			    "jne 0b                   \n\t"
			    "1:                       \n\t"
			    :"+r"(j)
			    :"r"(b), "r"(c), "g"(curr_words & (uint32)(~3))
			    :"%mm0", "%mm1", "%mm2", "%mm3", "memory", "cc");
#else
			for (j = 0; j < (curr_words & (uint32)(~3)); j += 4) {
				b[j+0] += c[j+0];
				b[j+1] += c[j+1];
				b[j+2] += c[j+2];
				b[j+3] += c[j+3];
			}
#endif
			for (; j < curr_words; j++)
				b[j] += c[j];

			block_offset += curr_words;
			contrib_offset = mp_modadd_1(contrib_offset,
						curr_words, contrib_words);
		}

		sp->contrib_array_offset = contrib_offset;
	}

#if defined(GCC_ASM32X) && defined(HAS_MMX)
	ASM_G volatile("emms");
#endif
}

/*-------------------------------------------------------------------------*/
static void
prepare_sieve_line(sieve_prime_t *primes, uint32 num_primes, 
			int64 xmin, int32 y)
{
	uint32 i, j, k;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint16 *contrib_array = curr_prime->contrib_array;
		uint32 contrib_array_size = curr_prime->contrib_array_size;

		memset(contrib_array, 0, contrib_array_size * sizeof(uint16));

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			uint32 num_roots = sp->num_roots;
			uint32 power = sp->power;
			uint16 contrib = sp->sieve_contrib;
			uint32 xmin_mod, y_mod;
			int32 tmpval;

			tmpval = xmin % (int32)power;
			if (tmpval < 0)
				tmpval += power;
			xmin_mod = tmpval;

			tmpval = y % (int32)power;
			if (tmpval < 0)
				tmpval += power;
			y_mod = tmpval;

			for (k = 0; k < num_roots; k++) {
				sieve_root_t *r = sp->roots + k;
				uint32 start = r->start;
				uint32 step = r->step;
				uint32 resclass = r->resclass;

				start = mp_modsub_1(start, xmin_mod, power);
				start = mp_modsub_1(start, resclass * y_mod % 
						    power, power);

				if (step != power)
					start = start % step;

				while (start < contrib_array_size) {
					contrib_array[start] += contrib;
					start += step;
				}
			}
		}

		curr_prime->contrib_array_offset = 0;

		for (j = 1; j < UNROLL; j++) {
			memcpy(contrib_array + j * contrib_array_size, 
				contrib_array, 
				contrib_array_size * sizeof(uint16));
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
do_sieving_line(sieve_prime_t *primes, uint32 num_primes,
			int64 xmin, int64 xmax, int32 y, 
			root_heap_t *heap, uint16 *block)
{
	uint32 i;
	uint32 worst_score = 0;

	prepare_sieve_line(primes, num_primes, xmin, y);

	if (heap->num_entries == heap->max_entries)
		worst_score = heap->entries[0].score;

	while (xmin < xmax) {
		int64 curr_line = xmax - xmin;
		int64 block_size = curr_line;

		if (curr_line > DEFAULT_BLOCK_SIZE &&
		    curr_line < 3 * DEFAULT_BLOCK_SIZE / 2)
			block_size = block_size / 2;

		block_size += UNROLL - block_size % UNROLL;
		block_size = MIN(block_size, DEFAULT_BLOCK_SIZE);

		sieve_one_block(block, (uint32)block_size, 
				primes, num_primes);

		for (i = 0; i < (uint32)block_size; i++) {
			uint32 score = block[i];
			if (score >= worst_score) {
				int64 off = xmin + i;

				if (off >= xmax)
					break;

				save_rotation(heap, off, y, score);
				if (heap->num_entries == heap->max_entries)
					worst_score = heap->entries[0].score;
			}
		}

		xmin += block_size;
	}
}

/*-------------------------------------------------------------------------*/
static void
prepare_sieve_lattice(sieve_prime_t *primes, uint32 num_primes,
			int64 xmin, int32 y, uint32 lattice_size, 
			uint32 lattice_class, uint16 *inv)
{
	uint32 i, j, k, m;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint16 *contrib_array = curr_prime->contrib_array;
		uint32 contrib_array_size = curr_prime->contrib_array_size;

		memset(contrib_array, 0, contrib_array_size * sizeof(uint16));

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			uint32 power = sp->power;
			uint32 num_roots = sp->num_roots;
			uint16 contrib = sp->sieve_contrib;
			uint32 mod_lattice = lattice_size % power;
			uint32 default_step;
			uint32 xmin_mod, y_mod;
			int32 tmpval;

			if (mod_lattice == 0)
				continue;

			default_step = power / mp_gcd_1(mod_lattice, power);

			tmpval = xmin % (int32)power;
			if (tmpval < 0)
				tmpval += power;
			xmin_mod = tmpval;

			tmpval = y % (int32)power;
			if (tmpval < 0)
				tmpval += power;
			y_mod = tmpval;

			for (k = 0; k < power; k++)
				inv[k] = power;

			for (k = lattice_class % power, m = 0; 
						m < power; m++) {
				if (inv[k] == power)
					inv[k] = m;
				k = mp_modadd_1(k, mod_lattice, power);
			}

			for (k = 0; k < num_roots; k++) {
				sieve_root_t *r = sp->roots + k;
				uint32 start = r->start;
				uint32 step = r->step;
				uint32 resclass = r->resclass;

				start = mp_modsub_1(start, xmin_mod, power);
				start = mp_modsub_1(start, resclass * y_mod % 
						    power, power);

				if (step == power) {
					start = inv[start];
					if (start == power)
						continue;
					step = default_step;
				}
				else {
					uint32 mult = power / step;

					start = start % step;
					start = inv[mult * start];
					if (start == power)
						continue;
					step /= mp_gcd_1(mod_lattice, step);
					start = start % step;
				}

				while (start < contrib_array_size) {
					contrib_array[start] += contrib;
					start += step;
				}
			}
		}

		curr_prime->contrib_array_offset = 0;

		for (j = 1; j < UNROLL; j++) {
			memcpy(contrib_array + j * contrib_array_size, 
				contrib_array, 
				contrib_array_size * sizeof(uint16));
		}
	}
}

/*-------------------------------------------------------------------------*/
static void
do_sieving_lattice(sieve_prime_t *primes, uint32 num_primes,
			int64 xmin, int64 xmax, int32 y0,
			uint32 lattice_size, root_heap_t *root_heap, 
			root_heap_t *lattice_heap, uint16 *block)
{
	uint32 i, j;
	uint32 worst_score = 0;
	int64 new_xmin = xmin / (int64)lattice_size - 1;
	int64 new_xmax = xmax / (int64)lattice_size + 1;

	if (root_heap->num_entries == root_heap->max_entries)
		worst_score = root_heap->entries[0].score;

	for (i = 0; i < lattice_heap->num_entries; i++) {

		uint32 resclass = lattice_heap->entries[i].x;
		int32 y = y0 + lattice_heap->entries[i].y;
		uint32 bias = lattice_heap->entries[i].score;
		int64 tmp_xmin = new_xmin;
		int64 tmp_xmax = new_xmax;

		prepare_sieve_lattice(primes, num_primes, 
					tmp_xmin * lattice_size, y, 
					lattice_size, resclass, block);

		while (tmp_xmin < tmp_xmax) {

			uint64 curr_line = tmp_xmax - tmp_xmin;
			uint64 block_size = curr_line;

			if (curr_line > DEFAULT_BLOCK_SIZE &&
			    curr_line < 3 * DEFAULT_BLOCK_SIZE / 2)
				block_size = block_size / 2;

			block_size += UNROLL - block_size % UNROLL;
			block_size = MIN(block_size, 
					DEFAULT_BLOCK_SIZE);

			sieve_one_block(block, (uint32)block_size, 
					primes, num_primes);

			for (j = 0; j < (uint32)block_size; j++) {
				uint32 score = bias + block[j];
				if (score >= worst_score) {
					int64 off = resclass + (tmp_xmin +
							j) * lattice_size;

					if (off < xmin)
						continue;

					if (off >= xmax)
						break;

					save_rotation(root_heap, off, y, score);

					if (root_heap->num_entries ==
						root_heap->max_entries) {
						worst_score = 
							root_heap->entries[
							      0].score;
					}
				}
			}

			tmp_xmin += block_size;
		}
	}
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_size(int64 line_size)
{
	if (line_size < 10000000)
		return 2*2*2*2*3*3*5;
	else if (line_size < 50000000)
		return 2*2*2*2*2*3*3*5;
	else if (line_size < 200000000)
		return 2*2*2*2*2*3*3*3*5;
	else if (line_size < 500000000)
		return 2*2*2*2*2*2*3*3*3*5;
	else if (line_size < (int64)1 * 1000000000)
		return 2*2*2*2*2*2*3*3*3*5*7;
	else if (line_size < (int64)3 * 1000000000)
		return 2*2*2*2*2*2*3*3*3*5*5*7;
	else if (line_size < (int64)10 * 1000000000)
		return 2*2*2*2*2*2*2*3*3*3*5*5*7;
	else if (line_size < (int64)30 * 1000000000)
		return 2*2*2*2*2*2*2*3*3*3*3*5*5*7;
	else if (line_size < (int64)200 * 1000000000)
		return 2*2*2*2*2*2*2*3*3*3*3*5*5*7*7;
	else if (line_size < (int64)1000 * 1000000000)
		return 2*2*2*2*2*2*2*3*3*3*3*5*5*5*7*7;
	else if (line_size < (int64)5000 * 1000000000)
		return 2*2*2*2*2*2*2*3*3*3*3*5*5*5*7*7*7;

	return 2*2*2*2*2*2*2*2*3*3*3*3*5*5*5*7*7*7;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			uint32 lattice_size, sieve_prime_t *lattice_primes)
{
	uint32 i, j;
	uint32 num_lattice_primes = 0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 num_lattice_powers = 0;

		if (lattice_size % curr_prime->prime)
			break;

		for (j = 0; j < num_powers; j++, num_lattice_powers++) {
			sieve_power_t *sp = curr_prime->powers + j;
			if (lattice_size % sp->power)
				break;
		}

		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 
						num_lattice_powers;
		num_lattice_primes++;
	}

	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
static void
find_lattices_core(sieve_prime_t *primes, uint32 num_primes,
			int32 y, uint32 lattice_size,
			root_heap_t *lattice_heap, uint16 *block)
{
	sieve_prime_t lattice_primes[10];
	uint32 num_lattice_primes;

	num_lattice_primes = find_lattice_primes(primes, num_primes, 
					lattice_size, lattice_primes);

	do_sieving_line(lattice_primes, num_lattice_primes,
			(int64)0, (int64)lattice_size, y,
			lattice_heap, block);
}

/*-------------------------------------------------------------------------*/
static void
find_lattices(root_sieve_t *rs, int32 y, uint32 lattice_size)
{
	sieve_prime_t lattice_primes[10];
	root_heap_t *lattice_heap = &rs->lattice_heap;
	root_heap_t *tmp_heap = &rs->tmp_lattice_heap;
	uint16 *block = rs->sieve_block;
	uint32 tmp_lattice_size;
	uint32 num_lattice_primes;

	lattice_heap->num_entries = 0;
	lattice_heap->max_entries = LATTICE_HEAP_SIZE;
	if (lattice_size < 1000)
		lattice_heap->max_entries /= 2;

	if (lattice_size < 200000) {
		find_lattices_core(rs->primes, rs->num_primes, y,
				   lattice_size, lattice_heap, block);
		return;
	}

	num_lattice_primes = find_lattice_primes(rs->primes, 
					rs->num_primes, lattice_size, 
					lattice_primes);

	tmp_lattice_size = find_lattice_size((int64)lattice_size);

	tmp_heap->max_entries = LATTICE_HEAP_SIZE;
	tmp_heap->num_entries = 0;
	find_lattices_core(lattice_primes, num_lattice_primes, y,
				tmp_lattice_size, tmp_heap, block);

	do_sieving_lattice(lattice_primes, num_lattice_primes,
			   (int64)0, (int64)lattice_size, 0,
			   tmp_lattice_size, lattice_heap, 
			   tmp_heap, block);
}

/*-------------------------------------------------------------------------*/
static void
root_sieve_x(root_sieve_t *rs, int64 xmin, int64 xmax, int32 y)
{
	int64 line_size = xmax - xmin;

	if (line_size < 10000000) {
		rs->root_heap.cutoffs[0].x = 50000;
		rs->root_heap.cutoffs[1].x = 1000000;
	}
	else if (line_size < 200000000) {
		rs->root_heap.cutoffs[0].x = 1000000;
		rs->root_heap.cutoffs[1].x = 5000000;
	}
	else {
		rs->root_heap.cutoffs[0].x = 10000000;
		rs->root_heap.cutoffs[1].x = 50000000;
	}

	if (line_size <= 1000000) {
		do_sieving_line(rs->primes, rs->num_primes, 
				xmin, xmax, y, &rs->root_heap,
				rs->sieve_block);
	}
	else {
		uint32 lattice_size = find_lattice_size(line_size);

		if (abs(y) < 50) {
			int64 small_xmin = MAX(-250000, xmin);
			int64 small_xmax = MIN(250000, xmax);

			if (small_xmin < 0 && small_xmax > 0) {
				do_sieving_line(rs->primes, rs->num_primes, 
						small_xmin, small_xmax, y, 
						&rs->root_heap, 
						rs->sieve_block);
			}
		}

		find_lattices(rs, y, lattice_size);

		do_sieving_lattice(rs->primes, rs->num_primes, 
				xmin, xmax, 0, lattice_size, 
				&rs->root_heap, &rs->lattice_heap, 
				rs->sieve_block);
	}
}

/*-------------------------------------------------------------------------*/
static uint32
root_sieve_xy(root_sieve_t *rs, double max_norm,
		double *a, uint32 deg, double p, double d,
		int32 y0, int32 y_inc)
{
	double a1 = a[1];
        double a2 = a[2];
	int64 xmin = -10000;
	int64 xmax = 10000;
	uint32 lines_done = 0;

	while (1) {
		a[1] = a1 - y0 * d;
		a[2] = a2 + y0 * p;
		compute_line_size(max_norm, a, deg, p, d, 
				  xmin, xmax, &xmin, &xmax);

//		printf("%d %.0lf %.0lf\n", y0, (double)xmin, (double)xmax);
		if (xmax - xmin < 1000)
			break;

		root_sieve_x(rs, xmin, xmax, y0);
		y0 += y_inc;
		lines_done++;
	}

	a[1] = a1;
	a[2] = a2;
	return lines_done;
}

/*-------------------------------------------------------------------------*/
static void process_rotations(poly_stage2_t *data, 
			root_sieve_t *rs) {

	uint32 i;

	for (i = 0; i < rs->root_heap.num_entries; i++) {
		rotation_t *r = rs->root_heap.entries + i;

		if (rs->random_root_score - (rs->sieve_bias + 
			(double)r->score / LOG_SCALE_FACTOR) <
					ROOT_SCORE_COARSE_MIN) {
			optimize_final(r->x, r->y, data);
		}
	}
}

/*-------------------------------------------------------------------------*/
static uint32
root_sieve_run_core(poly_stage2_t *data, double alpha_proj)
{
	uint32 i;
	double dbl_p;
	double dbl_d;
	double dbl_sv[MAX_POLY_DEGREE + 1];
	double max_norm;
	uint32 deg = data->degree;
	uint32 lines_done = 0;

	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	root_sieve_t *rs = &s->root_sieve;

	dbl_p = mpz_get_d(c->gmp_p);
	dbl_d = mpz_get_d(c->gmp_d);
	for (i = 0; i <= deg; i++)
		dbl_sv[i] = mpz_get_d(c->gmp_a[i]);

	max_norm = data->max_norm * exp(-alpha_proj);
	init_sieve(c, rs, deg, -alpha_proj);

	rs->root_heap.cutoffs[0].y = 100;
	rs->root_heap.cutoffs[0].score = LOG_SCALE_FACTOR * 
					(rs->random_root_score - 
					 rs->sieve_bias - good_alpha[deg-4][0]);
	rs->root_heap.cutoffs[1].y = 1000;
	rs->root_heap.cutoffs[1].score = LOG_SCALE_FACTOR * 
					(rs->random_root_score - 
					 rs->sieve_bias - good_alpha[deg-4][1]);
	rs->root_heap.extra = data;
	rs->root_heap.default_cutoff = LOG_SCALE_FACTOR * 
					(rs->random_root_score - 
					 rs->sieve_bias - good_alpha[deg-4][2]);

	rs->root_heap.num_entries = 0;
	lines_done += root_sieve_xy(rs, max_norm, dbl_sv, 
					deg, dbl_p, dbl_d, 0, 1);
	process_rotations(data, rs);

	rs->root_heap.num_entries = 0;
	lines_done += root_sieve_xy(rs, max_norm, dbl_sv, 
					deg, dbl_p, dbl_d, -1, -1);
	process_rotations(data, rs);

	return lines_done;
}

/*-------------------------------------------------------------------------*/
void
root_sieve_run(poly_stage2_t *data, double alpha_proj)
{
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	curr_poly_t *c = &s->curr_poly;
	uint32 deg = data->degree;
	uint32 z;

	if (deg != 6) {
		root_sieve_run_core(data, alpha_proj);
		return;
	}

	z = 0;
	while (1) {
		if (root_sieve_run_core(data, alpha_proj) == 0)
			break;

		z++;
		mpz_add(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
		mpz_sub(c->gmp_a[2], c->gmp_a[2], c->gmp_d);
	}
	mpz_submul_ui(c->gmp_a[3], c->gmp_p, (mp_limb_t)(z + 1));
	mpz_addmul_ui(c->gmp_a[2], c->gmp_d, (mp_limb_t)(z + 1));

	z = 1;
	while (1) {
		if (root_sieve_run_core(data, alpha_proj) == 0)
			break;

		z++;
		mpz_sub(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
		mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_d);
	}
	mpz_addmul_ui(c->gmp_a[3], c->gmp_p, (mp_limb_t)z);
	mpz_submul_ui(c->gmp_a[2], c->gmp_d, (mp_limb_t)z);
}

/*-------------------------------------------------------------------------*/
void root_sieve_init(root_sieve_t *rs)
{
	uint32 i, j;
	uint32 p;
	uint32 num_primes;

	memset(rs, 0, sizeof(root_sieve_t));

	rs->root_heap.max_entries = ROOT_HEAP_SIZE;
	rs->root_heap.entries = (rotation_t *)xmalloc(ROOT_HEAP_SIZE * 
							sizeof(rotation_t));
	rs->lattice_heap.max_entries = LATTICE_HEAP_SIZE;
	rs->lattice_heap.entries = (rotation_t *)xmalloc(LATTICE_HEAP_SIZE * 
							sizeof(rotation_t));
	rs->tmp_lattice_heap.max_entries = LATTICE_HEAP_SIZE;
	rs->tmp_lattice_heap.entries = (rotation_t *)xmalloc(LATTICE_HEAP_SIZE *
							sizeof(rotation_t));

	rs->sieve_block = (uint16 *)aligned_malloc(MAX(DEFAULT_BLOCK_SIZE,
						      MAX_SIEVE_PRIME_POWER) *
						   sizeof(uint16), 64);

	for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
		p += prime_delta[i];
		if (p > MAX_SIEVE_PRIME)
			break;

		rs->random_root_score += log((double)p) / (p - 1);
	}

	num_primes = rs->num_primes = i;
	rs->primes = (sieve_prime_t *)xmalloc(num_primes * 
						sizeof(sieve_prime_t));

	for (i = p = 0; i < num_primes; i++) {

		uint32 num_powers;
		uint32 power;
		double dlog;
		sieve_prime_t *curr_prime = rs->primes + i;

		p += prime_delta[i];
		num_powers = 1;
		power = p;
		while (power * p < MAX_SIEVE_PRIME_POWER) {
			num_powers++;
			power *= p;
		}

		dlog = log((double)p) * p / (p + 1);
		curr_prime->prime = p;
		curr_prime->num_powers = num_powers;
		curr_prime->powers = (sieve_power_t *)xmalloc(num_powers *
						sizeof(sieve_power_t));
		curr_prime->contrib_array_size = power;
		curr_prime->contrib_array = (uint16 *)aligned_malloc(
						power * UNROLL * 
						sizeof(uint16), 64);

		for (j = 0, power = p; j < num_powers; j++, power *= p) {

			sieve_power_t *curr_power = curr_prime->powers + j;

			curr_power->power = power;
			curr_power->root_contrib = dlog / power;
			curr_power->sieve_contrib = (uint16)(0.5 +
						LOG_SCALE_FACTOR * 
						curr_power->root_contrib);
			curr_power->roots = (sieve_root_t *)xmalloc(power *
							sizeof(sieve_root_t));
		}
	}
}

/*-------------------------------------------------------------------------*/
void root_sieve_free(root_sieve_t *rs)
{
	uint32 i, j;

	for (i = 0; i < rs->num_primes; i++) {
		sieve_prime_t *prime = rs->primes + i;

		for (j = 0; j < prime->num_powers; j++) {
			free(prime->powers[j].roots);
		}
		aligned_free(prime->contrib_array);
	}
	free(rs->primes);
	free(rs->root_heap.entries);
	free(rs->lattice_heap.entries);
	free(rs->tmp_lattice_heap.entries);
	aligned_free(rs->sieve_block);
	memset(rs, 0, sizeof(root_sieve_t));
}
