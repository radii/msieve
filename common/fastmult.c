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

/* Fast Hartley Transform based convolution routines

   See K. Jayasuriya, "Multiprecision Arithmetic Using
   Fast Hartley Transforms", Appl. Math. and Comp. 75:239-251 (1996)

   This paper explains more about the implementation aspects
   of FHTs than any other resource I've seen
*/

#include <fastmult.h>

#ifndef M_PI
#define M_PI  3.1415926535897932384626433832795029
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865475244008443621048490
#endif

/*------------------------------------------------------------------------*/
void fastmult_info_init(fastmult_info_t *info) {

	memset(info, 0, sizeof(fastmult_info_t));

	/* we use a tricky scheme for fast rounding of
	   floating point numbers, and have to choose a
	   constant that depends on the mantissa size of
	   the floating point registers. This is 53 bits
	   for all processors except x86 machines, where
	   it can be 53 or 64 bits depending on the machine,
	   compiler, OS, whether the code is built in debug 
	   or release mode, and possibly on the alignment of 
	   certain extrasolar planets. To avoid ambiguity, 
	   force the precision to always be 53 bits */

	if (!dd_precision_is_ieee()) {
		info->precision_changed = 1;
		info->old_precision = dd_set_precision_ieee();
	}
	info->round_constant[0] = 6755399441055744.0;
	info->round_constant[1] = 6755399441055744.0;
}

/*------------------------------------------------------------------------*/
void fastmult_info_free(fastmult_info_t *info) {
	
	uint32 i;
	uint32 num_twiddle = MIN(info->log2_runlength, HUGE_TWIDDLE_CUTOFF);
	uint32 num_huge_twiddle = info->log2_runlength - num_twiddle;

	for (i = 0; i <= num_twiddle; i++) {
		free(info->twiddle[i]);
	}
	for (i = 0; i <= num_huge_twiddle; i++) {
		free(info->huge_twiddle[i].small);
		free(info->huge_twiddle[i].large);
	}

	if (info->precision_changed) {
		dd_clear_precision(info->old_precision);
	}
	memset(info, 0, sizeof(fastmult_info_t));
}

/*------------------------------------------------------------------------*/
static double *make_twiddle(int32 size) {

	int32 i;
	double arg = (2.0 * M_PI) / size;
	double *w = (double *)xmalloc(2 * (size / 4 - 1) * sizeof(double));

	for (i = 0; i < size / 4 - 1; i++) {
		w[2 * i] = cos((i+1) * arg);
		w[2 * i + 1] = sin((i+1) * arg);
	}

	return w;
}

/*------------------------------------------------------------------------*/
static void make_huge_twiddle(huge_twiddle_t *table, int32 entry) {

	int32 i;
	int32 power = HUGE_TWIDDLE_CUTOFF + entry;
	double arg = (2.0 * M_PI) / (1 << power);
	int32 small_power = (power - 2) / 2;
	int32 large_power = power - 2 - small_power;

	table->large = (double *)xmalloc((2 << large_power) * sizeof(double));
	table->small = (double *)xmalloc((2 << small_power) * sizeof(double));

	for (i = 0; i < (1 << small_power); i++) {
		table->small[2 * i] = cos(i * arg);
		table->small[2 * i + 1] = sin(i * arg);
	}

	for (i = 0; i < (1 << large_power); i++) {
		table->large[2 * i] = cos((i << small_power) * arg);
		table->large[2 * i + 1] = sin((i << small_power) * arg);
	}
}

/*------------------------------------------------------------------------*/
static void create_twiddle(fastmult_info_t *info, int32 power) {

	int32 i;

	if (power <= (int32)(info->log2_runlength))
		return;

	if (power > MAX_FHT_POWER) {
		printf("error: transform size 1 << %d too large\n", power);
		exit(-1);
	}
	for (i = 4; i <= power; i++) {
		if (i < HUGE_TWIDDLE_CUTOFF) {
			if (info->twiddle[i] == NULL)
				info->twiddle[i] = make_twiddle(1 << i);
		}
		else {
			int32 j = i - HUGE_TWIDDLE_CUTOFF;
			huge_twiddle_t *h = info->huge_twiddle + j;
			if (h->large == NULL || h->small == NULL) {
				make_huge_twiddle(h, j);
			}
		}
	}
	info->log2_runlength = power;
}

/*------------------------------------------------------------------------*/
static void fht8_dif(double *x) {

	double t0, t1, t2, t3, t4, t5, t6, t7;
	double t8, t9, t10, t11, t12, t13, t14, t15;

	t0 = x[0];
	t2 = x[2];
	t4 = x[4];
	t6 = x[6];

	t8 = t0 + t4;
	t12 = t0 - t4;
	t10 = t6 + t2;
	t14 = t6 - t2;

	t1 = x[1];
	t3 = x[3];
	t5 = x[5];
	t7 = x[7];

	t9 = t1 + t5;
	t13 = t1 - t5;
	t11 = t7 + t3;
	t15 = t7 - t3;
	t6 = M_SQRT1_2 * (t15 + t13);
	t7 = M_SQRT1_2 * (t15 - t13);

	t0 = t8 + t10;
	t2 = t8 - t10;
	t1 = t11 + t9;
	t3 = t11 - t9;

	x[0] = t0 + t1;
	x[1] = t0 - t1;
	x[2] = t2 + t3;
	x[3] = t2 - t3;

	t0 = t12 + t14;
	t2 = t12 - t14;
	t1 = t7 + t6;
	t3 = t7 - t6;

	x[4] = t0 + t1;
	x[5] = t0 - t1;
	x[6] = t2 + t3;
	x[7] = t2 - t3;
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dif(double *x, double *w, int32 size) {

	double t0, t1, t2, t3, t4, t5;
	double c, s;
	double *x0, *x1;
	int32 i = size / 4;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t3 + t1;
	x1[-i] = t0 - t2;
	x1[0]  = t3 - t1;
	i--;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		c = w[0];
		s = w[1];

		x0[-i] = t0 + t2;
		x0[i]  = t1 + t3;
		t4 = t0 - t2;
		t5 = t3 - t1;
		x1[-i] = c * t4 + s * t5;
		x1[i] = c * t5 - s * t4;

		w += 2;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dif_huge(double *x, 
			huge_twiddle_t *h, int32 power) {

	int32 size = 1 << power;
	double t0, t1, t2, t3, t4, t5;
	double c, cs, cl, s, ss, sl;
	double *x0, *x1;
	int32 i = size / 4;
	int32 twiddle_arg;
	int32 shift = (power - 2) / 2;
	int32 mask = (1 << shift) - 1;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t3 + t1;
	x1[-i] = t0 - t2;
	x1[0]  = t3 - t1;
	i--;
	twiddle_arg = 1;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		cl = h->large[2*(twiddle_arg >> shift)];
		sl = h->large[2*(twiddle_arg >> shift) + 1];
		cs = h->small[2*(twiddle_arg & mask)];
		ss = h->small[2*(twiddle_arg & mask) + 1];
		c = cl * cs - sl * ss;
		s = cl * ss + cs * sl;

		x0[-i] = t0 + t2;
		x0[i]  = t1 + t3;
		t4 = t0 - t2;
		t5 = t3 - t1;
		x1[-i] = c * t4 + s * t5;
		x1[i] = c * t5 - s * t4;

		twiddle_arg++;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_forward(double *x, fastmult_info_t *info, int32 power) {

	int32 size;

	if (power == 3) {
		fht8_dif(x);
		return;
	}

	size = 1 << power;
	if (power < HUGE_TWIDDLE_CUTOFF) {
		fht_butterfly_dif(x, info->twiddle[power], size);
	}
	else {
		fht_butterfly_dif_huge(x, 
			info->huge_twiddle + (power - HUGE_TWIDDLE_CUTOFF), 
			power);
	}

	fht_forward(x, info, power - 1);
	fht_forward(x + size / 2, info, power - 1);
}

/*------------------------------------------------------------------------*/
static void fht8_dit(double *x) {

	double t0, t1, t2, t3, t4, t5, t6, t7;
	double t8, t9, t10, t11;

	t0 = x[0];
	t1 = x[1];
	t2 = x[2];
	t3 = x[3];

	t4 = t0 + t1;
	t5 = t0 - t1;
	t6 = t2 + t3;
	t7 = t2 - t3;

	t0 = t4 + t6;
	t2 = t4 - t6;
	t1 = t5 - t7;
	t3 = t5 + t7;

	t4 = x[4];
	t5 = x[5];
	t6 = x[6];
	t7 = x[7];

	t8 = t4 + t5;
	t9 = t4 - t5;
	t10 = t6 + t7;
	t11 = t6 - t7;

	t4 = t8 + t10;
	t6 = t8 - t10;
	t5 = t9 - t11;
	t7 = t9 + t11;

	x[0] = t0 + t4;
	x[4] = t0 - t4;
	x[2] = t2 - t6;
	x[6] = t2 + t6;

	t0 = M_SQRT1_2 * (t5 - t7);
	t2 = M_SQRT1_2 * (t5 + t7);

	x[1] = t1 + t0;
	x[5] = t1 - t0;
	x[3] = t3 - t2;
	x[7] = t3 + t2;
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dit(double *x, double *w, int32 size) {

	double t0, t1, t2, t3, t4, t5;
	double c, s;
	double *x0, *x1;
	int32 i = size / 4;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t1 - t3;
	x1[-i] = t0 - t2;
	x1[0]  = t1 + t3;
	i--;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		c = w[0];
		s = w[1];

		t4 = c * t2 - s * t3;
		t5 = s * t2 + c * t3;

		x0[-i] = t0 + t4;
		x0[i]  = t1 - t5;
		x1[-i] = t0 - t4;
		x1[i] = t1 + t5;

		w += 2;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dit_huge(double *x, 
			huge_twiddle_t *h, int32 power) {

	int32 size = 1 << power;
	double t0, t1, t2, t3, t4, t5;
	double c, cs, cl, s, ss, sl;
	double *x0, *x1;
	int32 i = size / 4;
	int32 twiddle_arg;
	int32 shift = (power - 2) / 2;
	int32 mask = (1 << shift) - 1;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t1 - t3;
	x1[-i] = t0 - t2;
	x1[0]  = t1 + t3;
	i--;
	twiddle_arg = 1;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		cl = h->large[2*(twiddle_arg >> shift)];
		sl = h->large[2*(twiddle_arg >> shift) + 1];
		cs = h->small[2*(twiddle_arg & mask)];
		ss = h->small[2*(twiddle_arg & mask) + 1];
		c = cl * cs - sl * ss;
		s = cl * ss + cs * sl;

		t4 = c * t2 - s * t3;
		t5 = s * t2 + c * t3;

		x0[-i] = t0 + t4;
		x0[i]  = t1 - t5;
		x1[-i] = t0 - t4;
		x1[i] = t1 + t5;

		twiddle_arg++;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_inverse(double *x, fastmult_info_t *info, int32 power) {

	int32 size;

	if (power == 3) {
		fht8_dit(x);
		return;
	}

	size = 1 << power;
	fht_inverse(x, info, power - 1);
	fht_inverse(x + size / 2, info, power - 1);

	if (power < HUGE_TWIDDLE_CUTOFF) {
		fht_butterfly_dit(x, info->twiddle[power], size);
	}
	else {
		fht_butterfly_dit_huge(x, 
			info->huge_twiddle + (power - HUGE_TWIDDLE_CUTOFF),
			power);
	}
}

/*------------------------------------------------------------------------*/
static void square_hermitian(double *x, int32 size) {

	int32 p;
	double r = 1.0 / (2 * size);

	x[0] = x[0] * x[0] * 2 * r;
	x[1] = x[1] * x[1] * 2 * r;

	for (p = 1; p < size/2; p = p * 2) {

		int32 i;
		int32 lo = 2 * p;
		int32 hi = 4 * p - 1;

		for (i = 0; i < p; i++, lo++, hi--) {
			double d1 = 2 * x[lo] * x[hi];
			double d2 = x[lo] * x[lo] - x[hi] * x[hi];
			x[lo] = (d1 + d2) * r;
			x[hi] = (d1 - d2) * r;
		}
	}
}

/*------------------------------------------------------------------------*/
static void mul_hermitian(double *x, double *y, int32 size) {

	int32 p;
	double r = 1.0 / (2 * size);

	x[0] = x[0] * y[0] * 2 * r;
	x[1] = x[1] * y[1] * 2 * r;

	for (p = 1; p < size/2; p = p * 2) {

		int32 i;
		int32 lo = 2 * p;
		int32 hi = 4 * p - 1;

		for (i = 0; i < p; i++, lo++, hi--) {
			double d1 = x[lo] * y[hi] + x[hi] * y[lo];
			double d2 = x[lo] * y[lo] - x[hi] * y[hi];
			x[lo] = (d1 + d2) * r;
			x[hi] = (d1 - d2) * r;
		}
	}
}

/*------------------------------------------------------------------------*/
static int32 get_log2_fhtsize(int32 nwords, int32 *log2_base) {

	int32 i;
	int32 runlength;
	int32 bits;

	/* begin by assuming base 2^16 for the FHT multiply */

	bits = 16;
	runlength = 1;
	while ((1 << runlength) < 2 * nwords)
		runlength++;

	/* if a base of 2^x for x = {17,18,19,20} would halve
	   the runlength, and the absolute maximum size of elements
	   of the product is not too much larger than the number
	   of bits in a double-precision mantissa, then use a base
	   larger than 2^16 */

	for (i = 17; i <= 20; i++) {
		if ( 2 * (i - 1) + (runlength - 1) <= 56 &&
		     (int32)(32.0 * nwords / i + 1) < 
		     			(1 << (runlength - 1))) {
			*log2_base = i;
			return runlength - 1;
		}
	}

	*log2_base = bits;
	return runlength;
}

/*------------------------------------------------------------------------*/
static void packed_to_double(uint32 *p, int32 pwords,
			 double *d, int32 fftsize,
			 int32 log2_base) {
	int32 i, j;

	uint64 accum = 0;		/* for operations on p */
	uint32 curr_bits = 0;
	uint32 shift = log2_base;
	uint32 mask = (1 << shift) - 1;

	int32 base = 1 << shift;	/* for operations on d */
	int32 half_base = base / 2;
	int32 word;
	int32 borrow;

	/* pull 32 bits at a time out of p and deposit 'shift'
	   bits at a time into d, converting to balanced form
	   along the way */

	i = j = borrow = 0;
	while (i < pwords) {

		if (curr_bits < shift) {
			accum |= (uint64)(p[i++]) << curr_bits;
			curr_bits += 32;
		}

		word = (int32)(accum & mask) + borrow;
		borrow = 0;
		accum >>= shift;
		curr_bits -= shift;

		if (word >= half_base) {
			word -= base;
			borrow = 1;
		}
		d[j++] = (double)word;
	}

	/* empty the accumulator */

	while (curr_bits) {
		word = (int32)(accum & mask) + borrow;
		borrow = 0;
		accum >>= shift;
		curr_bits -= MIN(curr_bits, shift);

		if (word >= half_base) {
			word -= base;
			borrow = 1;
		}
		d[j++] = (double)word;
	}
	d[j++] = borrow;

	/* pad with zeros */

	for (; j < fftsize; j++)
		d[j] = 0.0;
}

/*------------------------------------------------------------------------*/
static void double_to_packed(uint32 *p, int32 pwords,
				double *d, int32 log2_base,
				fastmult_info_t *info) {

	int32 i, j;

	double base = 1 << log2_base;	/* for operations on d */
	double recip_base = 1.0 / base;
	double carry = 0.0;
	double round0 = info->round_constant[0];
	double round1 = info->round_constant[1];
	double dword;

	uint64 accum = 0;		/* for operations on p */
	uint32 curr_bits = 0;
	uint32 shift = log2_base;

	/* pull 'shift' bits at a time out of d and deposit
	   32 bits at a time into p, converting from balanced
	   representation. Conversion stops when the expected
	   number of words in p have been found */

	i = j = 0;
	while (j < pwords) {
		dword = d[i++] + carry + round0 - round1;
		carry = dword * recip_base + round0 - round1;
		dword = dword - base * carry;
		if (dword < 0.0) {
			dword += base;
			carry--;
		}

		accum |= (uint64)((int32)dword) << curr_bits;
		curr_bits += shift;
		if (curr_bits >= 32) {
			p[j++] = (uint32)accum;
			accum >>= 32;
			curr_bits -= 32;
		}
	}

	/* check that the carry out from the product is zero.
	   The floating point array for holding the product
	   was allocated with one extra word, which should
	   cancel out the carry generated by the product being
	   in balanced form */

	dword = d[i] + carry + round0 - round1;
	carry = dword * recip_base + round0 - round1;
	dword = dword - base * carry;
	if (carry != 0 || dword != 0) {
		printf("error: overflow %lf carry %lf\n", dword, carry);
		exit(-1);
	}
}

/*------------------------------------------------------------------------*/
static void fht_square(int32 power, uint32 *a, int32 awords,
			uint32 *prod, int32 log2_base, 
			fastmult_info_t *info) {

	int32 runlength = 1 << power;
	double *factor1 = (double *)xmalloc(runlength * sizeof(double));

	packed_to_double(a, awords, factor1, runlength, log2_base);
	fht_forward(factor1, info, power);
	square_hermitian(factor1, runlength);
	fht_inverse(factor1, info, power);
	double_to_packed(prod, 2 * awords, factor1, log2_base, info);
	free(factor1);
}

/*------------------------------------------------------------------------*/
static void fht_mul(int32 power, uint32 *a, int32 awords,
		uint32 *b, int32 bwords, uint32 *prod,
		int32 log2_base, fastmult_info_t *info) {

	int32 runlength = 1 << power;
	double *factor1 = (double *)xmalloc(runlength * sizeof(double));
	double *factor2 = (double *)xmalloc(runlength * sizeof(double));

	packed_to_double(a, awords, factor1, runlength, log2_base);
	packed_to_double(b, bwords, factor2, runlength, log2_base);
	fht_forward(factor1, info, power);
	fht_forward(factor2, info, power);
	mul_hermitian(factor1, factor2, runlength);
	fht_inverse(factor1, info, power);
	double_to_packed(prod, awords + bwords, factor1, log2_base, info);
	free(factor1);
	free(factor2);
}

/*------------------------------------------------------------------------*/
void fastmult(uint32 *a, uint32 awords,
		uint32 *b, uint32 bwords,
		uint32 *prod, fastmult_info_t *info) {

	/* the product can never be more than awords+bwords
	   in size; however, we use balanced representation,
	   and it's possible for 'a' and 'b' to need one bit
	   more each. If they both need this extra bit, the
	   product will require an extra word. Failing to add
	   the extra word when awords+bwords is a power of 2
	   will mean that the top bit of the product wraps 
	   around to the bottom, corrupting the multiply */

	int32 log2_base;
	int32 power;

	power = get_log2_fhtsize((int32)(awords + bwords + 1),
				&log2_base);

	create_twiddle(info, power);
	if (a == b && awords == bwords) {
		fht_square(power, a, (int32)awords, prod, log2_base, info);
	}
	else {
		fht_mul(power, a, (int32)awords, b, (int32)bwords, 
				prod, log2_base, info);
	}
}
