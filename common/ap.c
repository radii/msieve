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

#include <mp_int.h>
#include <ap.h>

/*---------------------------------------------------------------*/
static void ap_check_nwords(ap_t *a, uint32 nwords) {

	if (a->num_alloc < nwords) {
		a->num_alloc = nwords + 100;
		a->val = (uint32 *)xrealloc(a->val, a->num_alloc * 
						sizeof(uint32));
	}
}

/*---------------------------------------------------------------*/
void ap_copy(ap_t *src, ap_t *dest) {

	if (src == dest)
		return;
	ap_check_nwords(dest, src->nwords);

	memcpy(dest->val, src->val, src->nwords * sizeof(uint32));
	dest->nwords = src->nwords;
	dest->sign = src->sign;
}

/*---------------------------------------------------------------*/
void ap_mp2ap(mp_t *src, uint32 sign, ap_t *dest) {

	if (mp_is_zero(src)) {
		dest->nwords = 0;
		dest->sign = POSITIVE;
		return;
	}

	ap_check_nwords(dest, src->nwords);
	memcpy(dest->val, src->val, src->nwords * sizeof(uint32));
	dest->nwords = src->nwords;
	dest->sign = sign;
}

/*---------------------------------------------------------------*/
void ap_si2ap(uint32 i, uint32 sign, ap_t *dest) {

	if (i == 0) {
		dest->nwords = 0;
		dest->sign = POSITIVE;
		return;
	}

	ap_check_nwords(dest, 1);
	dest->val[0] = i;
	dest->nwords = 1;
	dest->sign = sign;
}

/*---------------------------------------------------------------*/
uint32 ap_bits(ap_t *a) {

	uint32 i, bits, mask, top_word;

	if (ap_is_zero(a))
		return 0;

	i = a->nwords;
	bits = 32 * i;
	top_word = a->val[i - 1];

#if defined(GCC_ASM32X) || defined(GCC_ASM64X)

	ASM_G("bsrl %1, %0": "=r"(mask) : "rm"(top_word) : "cc");
	bits -= 31 - mask;
#else
	mask = 0x80000000;
	if ((top_word >> 16) == 0) {
		mask = 0x8000;
		bits -= 16;
	}

	while ( !(top_word & mask) ) {
		bits--;
		mask >>= 1;
	}
#endif

	return bits;
}

/*---------------------------------------------------------------*/
static void ap_add_abs(ap_t *a, ap_t *b, ap_t *sum) {

	/* a->nwords is assumed >= b->nwords */

	uint32 min_words, max_words;
	uint32 i;
	uint32 carry = 0;
	uint32 acc;

	max_words = a->nwords;
	min_words = b->nwords;
	ap_check_nwords(sum, max_words + 1);

	for (i = 0; i < min_words; i++) {
		acc = a->val[i] + carry;
		carry = (acc < a->val[i]);
		sum->val[i] = acc + b->val[i];
		carry += (sum->val[i] < acc);
	}

	for (; i < max_words; i++) {
		acc = a->val[i] + carry;
		carry = (acc < a->val[i]);
		sum->val[i] = acc;
	}

	if (carry)
		sum->val[i++] = carry;

	sum->nwords = num_nonzero_words(sum->val, i);
}

/*---------------------------------------------------------------*/
static void ap_sub_abs(ap_t *a, ap_t *b, ap_t *diff) {

	/* a->nwords is assumed >= b->nwords */

	uint32 min_words, max_words;
	uint32 i;
	uint32 borrow = 0;
	uint32 acc;

	max_words = a->nwords;
	min_words = b->nwords;
	ap_check_nwords(diff, max_words);

	for (i = 0; i < min_words; i++) {
		acc = a->val[i] - borrow;
		borrow = (acc > a->val[i]);
		diff->val[i] = acc - b->val[i];
		borrow += (diff->val[i] > acc);
	}

	for (; i < max_words; i++) {
		acc = a->val[i] - borrow;
		borrow = (acc > a->val[i]);
		diff->val[i] = acc;
	}

	diff->nwords = num_nonzero_words(diff->val, max_words);
}

/*---------------------------------------------------------------*/
void ap_add(ap_t *a, ap_t *b, ap_t *sum) {

	if (ap_is_zero(a)) {
		ap_copy(b, sum);
		return;
	}
	if (ap_is_zero(b)) {
		ap_copy(a, sum);
		return;
	}
	
	switch(2 * a->sign + b->sign) {

	case 2*POSITIVE + POSITIVE:
	case 2*NEGATIVE + NEGATIVE:
		if (ap_cmp_abs(a, b) >= 0)
			ap_add_abs(a, b, sum);
		else
			ap_add_abs(b, a, sum);
		sum->sign = a->sign;
		break;

	case 2*POSITIVE + NEGATIVE:
		if (ap_cmp_abs(a, b) >= 0) {
			ap_sub_abs(a, b, sum);
			sum->sign = POSITIVE;
		}
		else {
			ap_sub_abs(b, a, sum);
			sum->sign = NEGATIVE;
		}
		break;

	case 2*NEGATIVE + POSITIVE:
		if (ap_cmp_abs(a, b) > 0) {
			ap_sub_abs(a, b, sum);
			sum->sign = NEGATIVE;
		}
		else {
			ap_sub_abs(b, a, sum);
			sum->sign = POSITIVE;
		}
		break;
	}
}


/*---------------------------------------------------------------*/
void ap_sub(ap_t *a, ap_t *b, ap_t *diff) {

	if (ap_is_zero(a)) {
		ap_copy(b, diff);
		diff->sign = b->sign ^ 1;
		return;
	}
	if (ap_is_zero(b)) {
		ap_copy(a, diff);
		return;
	}
	
	switch(2 * a->sign + b->sign) {

	case 2*POSITIVE + POSITIVE:
		if (ap_cmp_abs(a, b) >= 0) {
			ap_sub_abs(a, b, diff);
			diff->sign = POSITIVE;
		}
		else {
			ap_sub_abs(b, a, diff);
			diff->sign = NEGATIVE;
		}
		break;

	case 2*NEGATIVE + NEGATIVE:
		if (ap_cmp_abs(a, b) > 0) {
			ap_sub_abs(a, b, diff);
			diff->sign = NEGATIVE;
		}
		else {
			ap_sub_abs(b, a, diff);
			diff->sign = POSITIVE;
		}
		break;

	case 2*POSITIVE + NEGATIVE:
	case 2*NEGATIVE + POSITIVE:
		if (ap_cmp_abs(a, b) >= 0)
			ap_add_abs(a, b, diff);
		else
			ap_add_abs(b, a, diff);
		diff->sign = a->sign;
		break;
	}
}

/*---------------------------------------------------------------*/
void ap_rshift(ap_t *a, uint32 shift, ap_t *res) {

	int32 i;
	int32 words = a->nwords;
	int32 start_word = shift / 32;
	uint32 word_shift = shift & 31;
	uint32 comp_word_shift = 32 - word_shift;

	if (start_word > words) {
		res->nwords = 0;
		res->sign = POSITIVE;
		return;
	}

	ap_check_nwords(res, (uint32)(words - start_word));

	if (word_shift == 0) {
		for (i = 0; i < (words-start_word); i++)
			res->val[i] = a->val[start_word+i];
	}
	else {
		for (i = 0; i < (words-start_word-1); i++) {
			res->val[i] = a->val[start_word+i] >> word_shift |
				a->val[start_word+i+1] << comp_word_shift;
		}
		res->val[i] = a->val[start_word+i] >> word_shift;
	}

	res->nwords = num_nonzero_words(res->val, (uint32)(words - start_word));
	res->sign = a->sign;
}

/*---------------------------------------------------------------*/
void ap_lshift(ap_t *a, uint32 shift, ap_t *res) {

	int32 i;
	uint32 words = a->nwords;
	uint32 start_word = shift / 32;
	uint32 word_shift = shift & 31;
	uint32 comp_word_shift = 32 - word_shift;

	if (ap_is_zero(a)) {
		res->nwords = 0;
		res->sign = POSITIVE;
		return;
	}

	ap_check_nwords(res, (uint32)(words + start_word + 1));

	if (word_shift == 0) {
		res->val[words + start_word] = 0;
		for (i = words - 1; (int32)i >= 0; i--)
			res->val[start_word + i] = a->val[i];
	}
	else {
		res->val[words + start_word] = 
					a->val[words - 1] >> comp_word_shift;
		for (i = words - 1; i; i--) {
			res->val[start_word + i] = a->val[i] << word_shift |
					a->val[i-1] >> comp_word_shift;
		}
		res->val[start_word + i] = a->val[i] << word_shift;
	}

	memset(res->val, 0, start_word * sizeof(uint32));
	res->nwords = num_nonzero_words(res->val, words + start_word + 1);
	res->sign = a->sign;
}

/*---------------------------------------------------------------*/
static void ap_addmul_1(uint32 *a, uint32 awords, uint32 b, uint32 *x) {

	uint32 carry = 0;

#if defined(GCC_ASM32A)

	uint32 tmp = awords;

	ASM_G(
	    "negl %0			\n\t"
	    "jz 1f			\n\t"
	    "0:				\n\t"
	    "movl (%2,%0,4), %%eax	\n\t"
	    "mull %4			\n\t"
	    "addl %1, %%eax		\n\t"
	    "adcl $0, %%edx		\n\t"
	    "addl %%eax, (%3,%0,4)	\n\t"
	    "movl %%edx, %1		\n\t"
	    "adcl $0, %1		\n\t"
	    "addl $1, %0		\n\t"
	    "jnz 0b			\n\t"
	    "1:				\n\t"
	    : "+r"(tmp), "+r"(carry)
	    : "r"(a + awords), "r"(x + awords), "m"(b)
	    : "%eax", "%edx", "cc", "memory");

#elif defined(MSC_ASM32A)

	ASM_M
	{
		push	ebx
		xor	ebx,ebx			; carry
		mov	ecx,awords		; negative loop count
		mov	esi,a			; pointer to source
		mov	edi,x			; pointer to destination
		lea	esi,[esi+ecx*4]
		lea	edi,[edi+ecx*4]
		neg	ecx
		jz	L1
	L0:	mov	eax,[esi+ecx*4]
		mul	b
		add	eax,ebx
		adc	edx,0
		add	[edi+ecx*4],eax
		mov	ebx,edx
		adc	ebx,0
		add	ecx,1
		jnz	L0
		mov	carry,ebx
	L1:	pop	ebx
	}

#elif defined(MSC_ASM64X)

	/* TODO: use 64-bit operations (but 
	   check if one 32-bit op is needed) */

	ASM_M		
	{				
		mov	r10,rcx		; entry rcx = *a, rdx = awoords
		mov	r11,r9		;       r8  =  b,  r9 = *x
		mov	rcx,rdx
		xor	r9,r9		; carry
		lea	r10,[r10+rcx*4]	; pointer to source
		lea	r11,[r11+rcx*4]	; pointer to destination
		neg	rcx		; note b is in r8 already
		jz	L1
	L0:	mov	eax,[r10+rcx*4]
		mul	r8d
		add	eax,r9d
		adc	edx,0
		add	[r11+rcx*4],eax
		mov	r9d,edx
		adc	r9d,0
		add	ecx,1
		jnz	L0
		mov	carry,r9d
	L1:	
	}

#else
	uint32 i;
	uint64 acc;

	for (i = 0; i < awords; i++) {
		acc = (uint64)a[i] * (uint64)b + 
		      (uint64)carry +
		      (uint64)x[i];
		x[i] = (uint32)acc;
		carry = (uint32)(acc >> 32);
	}
#endif
	
	x[awords] = carry;
}

/*---------------------------------------------------------------*/
static void ap_addmul(uint32 *a, uint32 awords,
		    uint32 *b, uint32 bwords,
		    uint32 *prod) {

	/* awords assumed >= bwords */

	uint32 i;

	for (i = 0; i < bwords; i++)
		ap_addmul_1(a, awords, b[i], prod + i);
}

/*---------------------------------------------------------------*/
void ap_mul(ap_t *a, ap_t *b, ap_t *prod, fastmult_info_t *info) {

	ap_t *c, *d;
	uint32 cwords, dwords, prod_words;

	if (ap_is_zero(a) || ap_is_zero(b)) {
		prod->nwords = 0;
		prod->sign = POSITIVE;
		return;
	}

	if (a->nwords > b->nwords) {
		c = a; d = b;
	}
	else {
		c = b; d = a;
	}
	cwords = c->nwords;
	dwords = d->nwords;
	prod_words = cwords + dwords;
	ap_check_nwords(prod, prod_words);

	if (cwords <= FFT_MIN_WORDS) {
		uint32 tmp[2 * FFT_MIN_WORDS];

		memset(tmp, 0, prod_words * sizeof(uint32));
		ap_addmul(c->val, cwords, d->val, dwords, tmp);
		memcpy(prod->val, tmp, prod_words * sizeof(uint32));
	}
	else if (dwords <= FFT_MIN_WORDS) {
		uint32 i, j;
		uint32 tmp[2 * FFT_MIN_WORDS] = {0};
		uint32 tmp_d_val[FFT_MIN_WORDS];
		ap_t tmp_d;
		uint32 mul_words = MIN(cwords, 2 * FFT_MIN_WORDS - dwords);

		if (prod == d) {
			tmp_d = *d;
			memcpy(tmp_d_val, d->val, dwords * sizeof(uint32));
			d = &tmp_d;
			d->val = tmp_d_val;
		}

		for (i = 0; i < cwords - mul_words; i += mul_words) {
			ap_addmul(c->val + i, mul_words, d->val, dwords, tmp);
			memcpy(prod->val + i, tmp, mul_words * sizeof(uint32));
			for (j = 0; j < dwords; j++)
				tmp[j] = tmp[j + mul_words];
			for (; j < 2 * FFT_MIN_WORDS; j++)
				tmp[j] = 0;
		}
		if (cwords - i < dwords)
			ap_addmul(d->val, dwords, c->val + i, cwords - i, tmp);
		else
			ap_addmul(c->val + i, cwords - i, d->val, dwords, tmp);

		memcpy(prod->val + i, tmp, 
				(cwords - i + dwords) * sizeof(uint32));
	}
	else {
		fastmult(c->val, cwords, d->val, dwords, prod->val, info);
	}

	prod->nwords = num_nonzero_words(prod->val, prod_words);
	prod->sign = c->sign ^ d->sign;
}

/*---------------------------------------------------------------*/
static void ap_mod_1(ap_t *num, ap_t *den, ap_t *res) {

	uint32 nwords = num->nwords;
	uint32 dwords = den->nwords;
	big_mp_t n;
	mp_t d, q, r;

	ap_check_nwords(res, dwords);
	if (dwords == 1) {
		res->val[0] = mp_mod_1_core(num->val, nwords, den->val[0]);
		res->nwords = (res->val[0]) ? 1 : 0;
		res->sign = num->sign;
		return;
	}

	mp_clear(&r);
	mp_clear(&d);
	d.nwords = dwords;
	memcpy(d.val, den->val, dwords * sizeof(uint32));

	while (nwords > 0) {
		uint32 i;
		uint32 chunk = MIN(nwords, MAX_MP_WORDS);
		
		for (i = 0; i < chunk; i++)
			n.val[i] = num->val[nwords - chunk + i];
		for (i = 0; i < r.nwords; i++)
			n.val[chunk + i] = r.val[i];
		for (i = chunk + r.nwords; i < 2 * MAX_MP_WORDS; i++)
			n.val[i] = 0;
		n.nwords = chunk + r.nwords;

		mp_divrem_core(&n, &d, &q, &r);
		nwords -= chunk;
	}

	ap_check_nwords(res, dwords);
	memcpy(res->val, r.val, dwords * sizeof(uint32));
	res->sign = num->sign;
	res->nwords = num_nonzero_words(res->val, dwords);
}

/*---------------------------------------------------------------*/
#define NUM_GUARD_BITS 32

void ap_recip(ap_t *a, ap_t *res, uint32 div_bits, fastmult_info_t *info) {

	uint32 curr_bits, new_bits;
	ap_t r2, a2;
	big_mp_t n;
	mp_t d, init_q, init_r;
	uint32 abits = ap_bits(a);
	uint32 prod_bits = abits + div_bits;

	/* this is a heavily modified version of the generalized
	   reciprocal algorithm from Crandall and Pomerance. In
           particular:

             - the precision is controlled adaptively, so that
	       the entire reciprocal process has an asymptotic 
	       latency of <= 4 full-precision multiplies

             - the iteration process can produce a reciprocal 
	       large enough to divide numbers with up to 
	       (div_bits+bits(a)) bits in one step. C&P specialize
	       to the case case of prod_bits = 2*bits(a)
	*/

	memset(&n, 0, sizeof(big_mp_t));
	mp_clear(&d);

	/* to get the initial approximation for use in the
	   Newton step, calculate 2^x / (a >> y), where the
	   numerator and denominator are chosen to be small 
	   enough for the quotient to be computed directly */

	curr_bits = MIN(abits, 32 * (MAX_MP_WORDS - 1));
	ap_rshift(a, abits - curr_bits, res);
	d.nwords = res->nwords;
	memcpy(d.val, res->val, d.nwords * sizeof(uint32));

	new_bits = MIN(prod_bits, curr_bits + 32 * (MAX_MP_WORDS - 1));
	n.nwords = new_bits / 32 + 1;
	n.val[n.nwords - 1] = 1 << (new_bits % 32);

	/* if x is large enough and y is zero, we have the
	   answer already */

	mp_divrem_core(&n, &d, &init_q, &init_r);
	ap_mp2ap(&init_q, POSITIVE, res);
	if (new_bits == prod_bits && curr_bits == abits)
		return;

	/* iterate until log2(answer) == div_bits */

	ap_init(&r2);
	ap_init(&a2);

	while (1) {
		/* each iteration will double the number of correct
		   bits in res. If doubling the precision will produce
		   more than div_bits correct bits, then reduce the
		   precision of the answer until doubling will provide
		   slightly more correct bits than we need */

		curr_bits = ap_bits(res);
		if (div_bits < 2 * curr_bits - NUM_GUARD_BITS) {
			ap_rshift(res, curr_bits - 
					(div_bits + NUM_GUARD_BITS) / 2, res);
			curr_bits = ap_bits(res);
		}

		/* square the previous answer. The number of bits
		   in the product will be the new precision level */

		ap_mul(res, res, &r2, info);
		new_bits = ap_bits(&r2);

		/* we have to get (a*res^2 >> new_bits) and 
		   (2 * previous_answer) to the current precision level. 
		   The latter is easy, and just needs a left shift. 
		   The former needs only the high-order bits of 'a' 
		   if the precision is low, or all of 'a' if it is high. */

		ap_lshift(res, new_bits - curr_bits + 1, res);
		if (abits <= new_bits) {
			ap_mul(&r2, a, &r2, info);
		}
		else {
			ap_rshift(a, abits - new_bits, &a2);
			ap_mul(&r2, &a2, &r2, info);
		}
		ap_rshift(&r2, ap_bits(&r2) - new_bits, &r2);

		/* compute the next approximation, and if it has
		   enough bits then we're done */

		ap_sub(res, &r2, res);
		new_bits = ap_bits(res);
		if (div_bits < new_bits) {
			ap_rshift(res, new_bits - div_bits - 1, res);
			break;
		}
	}

	ap_clear(&r2);
	ap_clear(&a2);
}

/*---------------------------------------------------------------*/
void ap_mod(ap_t *num, ap_t *den, ap_t *recip,
		ap_t *res, fastmult_info_t *info) {

	/* the algorithm is from Crandall and Pomerance, who
	   cite the Handbook of Applied Cryptography. Note that
	   we generalize the algorithm from C&P so that 
	   
	   - recip can be *any* size; in particular num can 
	     exceed den*den in size

	   - the division takes place log2(recip) bits at a time.
	     This lets calling code deal with huge operands in 
	     chunks of more manageable size 
	 */

	uint32 dbits;
	uint32 nbits1, nbits2;
	uint32 rbits;
	ap_t tmp;
	ap_t *curr_num;

	if (ap_is_zero(num) || ap_cmp_abs(num, den) == 0) {
		res->nwords = 0;
		res->sign = POSITIVE;
		return;
	}
	if (ap_cmp_abs(num, den) < 0) {
		ap_copy(num, res);
		return;
	}
	if (den->nwords <= MAX_MP_WORDS) {
		ap_mod_1(num, den, res);
		return;
	}

	/* perform as many reduction steps as are needed
	   to get the remainder to the neighborhood of
	   the correct result */

	ap_init(&tmp);
	dbits = ap_bits(den);
	rbits = ap_bits(recip);
	curr_num = num;
	nbits1 = ap_bits(curr_num);
	do {
		/* each iteration removes at most rbits bits from the
		   numerator. First compute an approximation to
		   the quotient, using MIN(rbits, bits(curr_num)) bits 
		   of recip and curr_num */

		if (nbits1 > rbits) {
			ap_rshift(curr_num, nbits1 - rbits, &tmp);
			ap_mul(&tmp, recip, &tmp, info);
		}
		else {
			ap_rshift(recip, rbits - nbits1, &tmp);
			ap_mul(curr_num, &tmp, &tmp, info);
		}

		/* compute the high-order bits of the quotient
		   and multiply by den */

		nbits2 = ap_bits(&tmp);
		if (nbits2 > nbits1 - dbits) {
			ap_rshift(&tmp, nbits2 - (nbits1 - dbits), &tmp);
		}
		ap_mul(&tmp, den, &tmp, info);

		/* compute num - quotient * den. Equalize the 
		   precision before subtracting */

		nbits2 = ap_bits(&tmp);
		if (nbits2 > nbits1)
			ap_rshift(&tmp, nbits2 - nbits1, &tmp);
		else
			ap_lshift(&tmp, nbits1 - nbits2, &tmp);
		ap_sub(curr_num, &tmp, res);
		curr_num = res;
		nbits1 = ap_bits(curr_num);
	} while (nbits1 > dbits + 1);

	/* compute the correct result from the approximation */

	while (ap_cmp_abs(res, den) >= 0) {
		if (res->sign == POSITIVE)
			ap_sub(res, den, res);
		else
			ap_add(res, den, res);
	}
	ap_clear(&tmp);
}
