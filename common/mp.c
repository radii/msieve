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

#include <mp.h>

/* silly hack: in order to expand macros and then
   turn them into strings, you need two levels of
   macro-izing */

#define _(x) #x
#define STRING(x) _(x)

/*---------------------------------------------------------------*/
static uint32 num_nonzero_words(uint32 *x, uint32 max_words) {

	/* return the index of the first nonzero word in
	   x, searching backwards from max_words */

	uint32 i;
	for (i = max_words; i && !x[i-1]; i--)
		;
	return i;
}

/*---------------------------------------------------------------*/
uint32 mp_bits(mp_t *a) {

	uint32 i, bits, mask, top_word;

	if (mp_is_zero(a))
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
double mp_log(mp_t *x) {

	uint32 i = x->nwords;

	switch(i) {
	case 0:
		return 0;
	case 1:
		return log((double)(x->val[0]));
	case 2:
		return log((double)(x->val[0]) + 
				MP_RADIX * x->val[1]);
	default:
		return 32 * (i-3) * M_LN2 + 
			log((double)(x->val[i-3]) + MP_RADIX * (
		     		((double)x->val[i-2] + MP_RADIX * 
				x->val[i-1])));
	}
}

/*---------------------------------------------------------------*/
void mp_add(mp_t *a, mp_t *b, mp_t *sum) {

	uint32 max_words;

#if defined(GCC_ASM32A)
	int32 index = -MAX_MP_WORDS;

	ASM_G volatile(
	    "xorl %%eax, %%eax					\n\t"
	    "0:                       				\n\t"
	    "movl " STRING(4*MAX_MP_WORDS) "(%1,%0,4), %%eax    \n\t"
	    "adcl " STRING(4*MAX_MP_WORDS) "(%2,%0,4), %%eax    \n\t"
	    "movl %%eax, " STRING(4*MAX_MP_WORDS) "(%3,%0,4)    \n\t"
	    "incl %0                  				\n\t"
	    "jnz 0b                   				\n\t"
	    :"+r"(index) 
	    :"r"(a->val), "r"(b->val), "r"(sum->val)
	    :"%eax", "memory", "cc");

	max_words = MAX(a->nwords, b->nwords);
	sum->nwords = max_words;
	if (max_words < MAX_MP_WORDS && sum->val[max_words])
		sum->nwords++;

#elif defined(GCC_ASM64X)
	int64 index = -MAX_MP_WORDS;

	ASM_G volatile(
	    "xorl %%eax, %%eax					\n\t"
	    "0:                       				\n\t"
	    "movl " STRING(4*MAX_MP_WORDS) "(%1,%0,4), %%eax    \n\t"
	    "adcl " STRING(4*MAX_MP_WORDS) "(%2,%0,4), %%eax    \n\t"
	    "movl %%eax, " STRING(4*MAX_MP_WORDS) "(%3,%0,4)    \n\t"
	    "incq %0                  				\n\t"
	    "jnz 0b                   				\n\t"
	    :"+r"(index) 
	    :"r"(a->val), "r"(b->val), "r"(sum->val)
	    :"%eax", "memory", "cc");

	max_words = MAX(a->nwords, b->nwords);
	sum->nwords = max_words;
	if (max_words < MAX_MP_WORDS && sum->val[max_words])
		sum->nwords++;

#elif defined(MSC_ASM32A)
	ASM_M
	{	
		xor	eax,eax
		mov	esi,a
		mov	edi,b
		mov	edx,sum
		lea	esi,[esi]a.val+4*MAX_MP_WORDS
		lea	edi,[edi]b.val+4*MAX_MP_WORDS
		lea	edx,[edx]sum.val+4*MAX_MP_WORDS
		mov	ecx,-MAX_MP_WORDS
	L0:	mov	eax,[esi+ecx*4]
		adc	eax,[edi+ecx*4]
		mov	[edx+ecx*4],eax
		inc	ecx
		jnz	L0
	}
	max_words = MAX(a->nwords, b->nwords);
	sum->nwords = max_words;
	if (max_words < MAX_MP_WORDS && sum->val[max_words])
		sum->nwords++;
#else
	uint32 i;
	uint32 carry = 0;
	uint32 acc;

	max_words = MAX(a->nwords, b->nwords);
	for (i = 0; i < max_words; i++) {
		acc = a->val[i] + carry;
		carry = (acc < a->val[i]);
		sum->val[i] = acc + b->val[i];
		carry += (sum->val[i] < acc);
	}
	if (carry)
		sum->val[i++] = carry;

	sum->nwords = i;
	for (; i < MAX_MP_WORDS; i++)
		sum->val[i] = 0;
#endif
}

/*---------------------------------------------------------------*/
void mp_add_1(mp_t *a, uint32 b, mp_t *sum) {

	uint32 max_words = a->nwords;
	uint32 i;
	uint32 carry = b;
	uint32 acc;

	for (i = 0; carry && i < max_words; i++) {
		acc = a->val[i] + carry;
		carry = (acc < a->val[i]);
		sum->val[i] = acc;
	}
	if (carry)
		sum->val[i++] = carry;

	for (; i < MAX_MP_WORDS; i++)
		sum->val[i] = a->val[i];

	sum->nwords = max_words;
	if (max_words < MAX_MP_WORDS && sum->val[max_words])
		sum->nwords++;
}

/*---------------------------------------------------------------*/
void signed_mp_add(signed_mp_t *a, signed_mp_t *b, signed_mp_t *sum) {

	switch(2 * a->sign + b->sign) {

	case 2*POSITIVE + POSITIVE:
	case 2*NEGATIVE + NEGATIVE:
		mp_add(&a->num, &b->num, &sum->num);
		sum->sign = a->sign;
		break;

	case 2*POSITIVE + NEGATIVE:
		if (mp_cmp(&a->num, &b->num) >= 0) {
			mp_sub(&a->num, &b->num, &sum->num);
			sum->sign = POSITIVE;
		}
		else {
			mp_sub(&b->num, &a->num, &sum->num);
			sum->sign = NEGATIVE;
		}
		break;

	case 2*NEGATIVE + POSITIVE:
		if (mp_cmp(&a->num, &b->num) > 0) {
			mp_sub(&a->num, &b->num, &sum->num);
			sum->sign = NEGATIVE;
		}
		else {
			mp_sub(&b->num, &a->num, &sum->num);
			sum->sign = POSITIVE;
		}
		break;
	}
}

/*---------------------------------------------------------------*/
void mp_sub(mp_t *a, mp_t *b, mp_t *diff) {

	uint32 max_words = a->nwords;

#if defined(GCC_ASM32A)

	int32 index = -MAX_MP_WORDS;
	ASM_G volatile(
	    "xorl %%eax, %%eax        \n\t"
	    "0:                       \n\t"
	    "movl " STRING(4*MAX_MP_WORDS) "(%1,%0,4), %%eax    \n\t"
	    "sbbl " STRING(4*MAX_MP_WORDS) "(%2,%0,4), %%eax    \n\t"
	    "movl %%eax, " STRING(4*MAX_MP_WORDS) "(%3,%0,4)    \n\t"
	    "incl %0                  \n\t"
	    "jnz 0b                   \n\t"
	    :"+r"(index)
	    :"r"(a->val), "r"(b->val), "r"(diff->val)
	    :"%eax", "memory", "cc");

#elif defined(GCC_ASM64X)

	int64 index = -MAX_MP_WORDS;
	ASM_G volatile(
	    "xorl %%eax, %%eax        \n\t"
	    "0:                       \n\t"
	    "movl " STRING(4*MAX_MP_WORDS) "(%1,%0,4), %%eax    \n\t"
	    "sbbl " STRING(4*MAX_MP_WORDS) "(%2,%0,4), %%eax    \n\t"
	    "movl %%eax, " STRING(4*MAX_MP_WORDS) "(%3,%0,4)    \n\t"
	    "incq %0                  \n\t"
	    "jnz 0b                   \n\t"
	    :"+r"(index)
	    :"r"(a->val), "r"(b->val), "r"(diff->val)
	    :"%eax", "memory", "cc");

#elif defined(MSC_ASM32A)
	ASM_M
	{
		xor	eax,eax
		mov	esi,a
		mov	edi,b
		mov	edx,diff
		lea	esi,[esi]a.val+4*MAX_MP_WORDS
		lea	edi,[edi]b.val+4*MAX_MP_WORDS
		lea	edx,[edx]diff.val+4*MAX_MP_WORDS
		mov	ecx,-MAX_MP_WORDS
	L0:	mov	eax,[esi+ecx*4]
		sbb	eax,[edi+ecx*4]
		mov	[edx+ecx*4],eax
		inc	ecx
		jnz	L0
	}
#else
	uint32 borrow = 0;
	uint32 acc;
	uint32 i;

	for (i = 0; i < max_words; i++) {
		acc = a->val[i] - borrow;
		borrow = (acc > a->val[i]);
		diff->val[i] = acc - b->val[i];
		borrow += (diff->val[i] > acc);
	}
	for (; i < MAX_MP_WORDS; i++)
		diff->val[i] = 0;
#endif

	diff->nwords = num_nonzero_words(diff->val, max_words);
}

/*---------------------------------------------------------------*/
void mp_sub_1(mp_t *a, uint32 b, mp_t *diff) {

	uint32 max_words = a->nwords;
	int32 i;
	uint32 borrow = b;
	uint32 acc;

	for (i = 0; borrow; i++) {
		acc = a->val[i] - borrow;
		borrow = (acc > a->val[i]);
		diff->val[i] = acc;
	}
	for (; i < MAX_MP_WORDS; i++)
		diff->val[i] = a->val[i];

	diff->nwords = num_nonzero_words(diff->val, max_words);
}

/*---------------------------------------------------------------*/
void signed_mp_sub(signed_mp_t *a, signed_mp_t *b, signed_mp_t *diff) {

	switch(2 * a->sign + b->sign) {

	case 2*POSITIVE + POSITIVE:
		if (mp_cmp(&a->num, &b->num) >= 0) {
			mp_sub(&a->num, &b->num, &diff->num);
			diff->sign = POSITIVE;
		}
		else {
			mp_sub(&b->num, &a->num, &diff->num);
			diff->sign = NEGATIVE;
		}
		break;

	case 2*NEGATIVE + NEGATIVE:
		if (mp_cmp(&a->num, &b->num) > 0) {
			mp_sub(&a->num, &b->num, &diff->num);
			diff->sign = NEGATIVE;
		}
		else {
			mp_sub(&b->num, &a->num, &diff->num);
			diff->sign = POSITIVE;
		}
		break;

	case 2*POSITIVE + NEGATIVE:
	case 2*NEGATIVE + POSITIVE:
		mp_add(&a->num, &b->num, &diff->num);
		diff->sign = a->sign;
		break;
	}
}

/*---------------------------------------------------------------*/
static void mp_addmul_1(mp_t *a, uint32 b, uint32 *x) {

	uint32 words = a->nwords;
	uint32 carry = 0;

#if defined(GCC_ASM32A)

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
	    : "+r"(words), "+r"(carry)
	    : "r"(a->val + words), "r"(x + words), "m"(b)
	    : "%eax", "%edx", "cc", "memory");

	words = a->nwords;

#elif defined(MSC_ASM32A)
	ASM_M
	{
		push	ebx
		xor	ebx,ebx
		mov	ecx,words
		mov	esi,a
		mov	edi,x
		lea	esi,[esi+4*ecx]a.val
		lea	edi,[edi+4*ecx]
		neg	ecx
		jz	L1
	L0:	mov	eax,[esi+ecx*4]
		mul	b
		add	eax,ebx
		adc	edx,0
		add	[edi+ecx*4],eax
		mov	ebx,edx
		adc	ebx,0
		add	ecx, 1
		jnz	L0
		mov	carry,ebx
	L1:	pop	ebx
	}
	words = a->nwords;

#else
	uint32 i;
	uint64 acc;

	for (i = 0; i < words; i++) {
		acc = (uint64)a->val[i] * (uint64)b + 
		      (uint64)carry +
		      (uint64)x[i];
		x[i] = (uint32)acc;
		carry = (uint32)(acc >> 32);
	}
#endif
	/* avoid writing the carry if zero. This is needed to
	   avoid a buffer overrun when callers are multiplying
	   integers whose product requires MP_MAX_WORDS words 
	   of storage. A side effect of this is that callers 
	   must guarantee that x[words] is initialized to zero 
	   when mp_addmul_1 is called */

	if (carry) {
		x[words] = carry;
	}
}

/*---------------------------------------------------------------*/
static uint32 mp_submul_1(uint32 *a, uint32 b, 
			uint32 words, uint32 *x) {

	uint32 carry = 0;

#if defined(GCC_ASM32A)
	ASM_G(
	    "negl %5			\n\t"
	    "jz 1f			\n\t"
	    "0:				\n\t"
	    "movl (%2,%5,4), %%eax	\n\t"
	    "mull %4			\n\t"
	    "addl %1, %%eax		\n\t"
	    "adcl $0, %%edx		\n\t"
	    "subl %%eax, (%3,%5,4)	\n\t"
	    "movl %%edx, %1		\n\t"
	    "adcl $0, %1		\n\t"
	    "addl $1, %0		\n\t"
	    "jnz 0b			\n\t"
	    "1:				\n\t"
	    : "+r"(words), "+r"(carry)
	    : "r"(a + words), "r"(x + words), "g"(b)
	    : "%eax", "%edx", "cc", "memory");

#elif defined(MSC_ASM32A)
	ASM_M
	{
		push	ebx
		xor	ebx,ebx
		mov	ecx,words
		mov	esi,a
		mov	edi,x
		lea	esi,[esi+4*ecx]
		lea	edi,[edi+4*ecx]
		neg	ecx
		jz	L1
	L0:	mov	eax,[esi+ecx*4]
		mul	b
		add	eax,ebx
		adc	edx,0
		sub	[edi+ecx*4],eax
		mov	ebx,edx
		adc	ebx,0
		add	ecx,1
		jnz	L0
		mov	carry,ebx
	L1:	pop	ebx
	}

#else
	uint32 i, j, k;

	for (i = 0; i < words; i++) {
		uint64 acc = (uint64)a[i] * (uint64)b + (uint64)carry;
		k = x[i];
		j = (uint32)acc;
		carry = (uint32)(acc >> 32);
		x[i] = k - j;
		carry += (x[i] > k);
	}
#endif
	return carry;
}

/*---------------------------------------------------------------*/
void mp_mul_1(mp_t *a, uint32 b, mp_t *x) {

	uint32 i;
	uint32 carry = 0;
	uint32 words = a->nwords;

	if (b == 0) {
		mp_clear(x);
		return;
	}

#if defined(GCC_ASM32A)
	ASM_G(
	    "negl %5			\n\t"
	    "jz 1f			\n\t"
	    "0:				\n\t"
	    "movl (%2,%5,4), %%eax	\n\t"
	    "mull %4			\n\t"
	    "addl %1, %%eax		\n\t"
	    "adcl $0, %%edx		\n\t"
	    "movl %%eax, (%3,%5,4)	\n\t"
	    "movl %%edx, %1		\n\t"
	    "addl $1, %5		\n\t"
	    "jnz 0b			\n\t"
	    "1:				\n\t"
	    : "+r"(words), "+r"(carry)
	    : "r"(a->val + words), "r"(x->val + words), "g"(b)
	    : "%eax", "%edx", "cc", "memory");

	 words = a->nwords;

#elif defined(MSC_ASM32A)
	ASM_M
	{
		push	ebx
		xor	ebx,ebx
		mov	ecx,words
		mov	esi,a
		mov	edi,x
		lea	esi,[esi+4*ecx]a.val
		lea	edi,[edi+4*ecx]x.val
		neg	ecx
		jz	L1
	L0:	mov	eax,[esi+ecx*4]
		mul	b
		add	eax,ebx
		adc	edx,0
		mov	[edi+ecx*4],eax
		mov	ebx,edx
		add	ecx,1
		jnz	L0
		mov	carry,ebx
	L1:	pop	ebx
	}
	words = a->nwords;

#else
	for (i = 0; i < words; i++) {
		uint64 acc = (uint64)a->val[i] * (uint64)b + (uint64)carry;
		x->val[i] = (uint32)acc;
		carry = (uint32)(acc >> 32);
	}
#endif

	if (carry) {
		x->val[words++] = carry;
	}
	x->nwords = words;
	for (i = words; i < MAX_MP_WORDS; i++)
		x->val[i] = 0;
}

/*---------------------------------------------------------------*/
void mp_mul(mp_t *a, mp_t *b, mp_t *prod) {

	uint32 i;
	mp_t *small = a;
	mp_t *large = b;

	if (small->nwords > large->nwords) {
		small = b;
		large = a;
	}
	if (small->nwords == 0) {
		mp_clear(prod);
		return;
	}

	mp_mul_1(large, small->val[0], prod);
	for (i = 1; i < small->nwords; i++) 
		mp_addmul_1(large, small->val[i], prod->val + i);
	
	prod->nwords = num_nonzero_words(prod->val, 
				MIN(a->nwords + b->nwords, MAX_MP_WORDS));
}

/*---------------------------------------------------------------*/
void signed_mp_mul(signed_mp_t *a, signed_mp_t *b, 
				signed_mp_t *prod) {

	mp_mul(&a->num, &b->num, &prod->num);
	prod->sign = a->sign ^ b->sign;
}

/*---------------------------------------------------------------*/
void mp_rshift(mp_t *a, uint32 shift, mp_t *res) {

	int32 i;
	int32 words = a->nwords;
	int32 start_word = shift / 32;
	uint32 word_shift = shift & 31;
	uint32 comp_word_shift = 32 - word_shift;

	if (start_word > words) {
		mp_clear(res);
		return;
	}

	if (word_shift == 0) {
		for (i = 0; i < (words-start_word); i++)
			res->val[i] = a->val[start_word+i];
	}
	else {
		for (i = 0; i < (words-start_word-1); i++)
			res->val[i] = a->val[start_word+i] >> word_shift |
				a->val[start_word+i+1] << comp_word_shift;
		res->val[i] = a->val[start_word+i] >> word_shift;
		i++;
	}

	for (; i < MAX_MP_WORDS; i++)
		res->val[i] = 0;

	res->nwords = num_nonzero_words(res->val, (uint32)(words - start_word));
}

/*---------------------------------------------------------------*/
uint32 mp_rjustify(mp_t *a, mp_t *res) {

	uint32 i, mask, shift, words;

	if (mp_is_zero(a))
		return 0;
	
	words = a->nwords;
	for (i = 0; i < words; i++) {
		if (a->val[i] != 0)
			break;
	}
	mask = a->val[i];
	shift = 32 * i;

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) 
	ASM_G("bsfl %1, %0": "=r"(i) : "rm"(mask));
#else
	for (i = 0; i < 32; i++) {
		if (mask & (1 << i))
			break;
	}
#endif
	shift += i;

	mp_rshift(a, shift, res);
	return shift;
}

/*---------------------------------------------------------------*/
void mp_divrem_core(big_mp_t *num, mp_t *denom, 
			mp_t *quot, mp_t *rem) {

	int32 i, j, k;
	uint32 shift, comp_shift;
	uint32 high_denom, low_denom;
	uint32 nacc[2*MAX_MP_WORDS+1];
	uint32 dacc[MAX_MP_WORDS];

	i = num->nwords;
	j = denom->nwords;
	high_denom = denom->val[j-1];

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) 
	ASM_G("bsrl %1, %0": "=r"(shift) : "rm"(high_denom) : "cc");
	shift = 31 - shift;
#else
	comp_shift = 0x80000000;
	shift = 0;
	if ((high_denom >> 16) == 0) {
		comp_shift = 0x8000;
		shift = 16;
	}

	while ( !(high_denom & comp_shift) ) {
		shift++;
		comp_shift >>= 1;
	}
#endif
	comp_shift = 32 - shift;
	if (shift) {
		for (k = j; k > 1; k--) {
			dacc[k-1] = denom->val[k-1] << shift |
					denom->val[k-2] >> comp_shift;
		}
		dacc[0] = denom->val[0] << shift;

		for (k = i; k > 1; k--) {
			nacc[k-1] = num->val[k-1] << shift |
					num->val[k-2] >> comp_shift;
		}
		nacc[0] = num->val[0] << shift;
		nacc[i] = num->val[i-1] >> comp_shift;
	}
	else {
		memcpy(dacc, denom->val, j * sizeof(uint32));
		memcpy(nacc, num->val, i * sizeof(uint32));
		nacc[i] = 0;
	}
	high_denom = dacc[j-1];
	low_denom = dacc[j-2];

	while (i + 1 > j) {
		uint32 q, r;
		uint32 high_num = nacc[i];
		uint32 low_num = nacc[i-1];
		uint32 check = nacc[i-2];

		if (high_num == high_denom) {

			/* the quotient will be 2^32, and must
			   always get corrected */

			q = 0xffffffff;
		}
		else {

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) 
			q = low_num;
			r = high_num;
			ASM_G("divl %2"
				: "+a"(q),"+d"(r)
				: "rm"(high_denom) : "cc");
#else
			uint64 acc = (uint64)high_num << 32 | (uint64)low_num;
			q = (uint32)(acc / high_denom);
			r = (uint32)(acc % high_denom);
#endif
			while ((uint64)low_denom * (uint64)q >
				((uint64)r << 32) + check) {
				q--;
				r += high_denom;
				if (r < high_denom)
					break;
			}
		}

		if (mp_submul_1(dacc, q, (uint32)j, nacc + i - j) 
							!= high_num) {
			uint32 carry = 0;
			for (q--, k = 0; k < j; k++) {
				uint32 sum = nacc[i-j+k] + carry;
				carry = (sum < nacc[i-j+k]);
				nacc[i-j+k] = sum + dacc[k];
				carry += (nacc[i-j+k] < sum);
			}
		}

		if (i - j < MAX_MP_WORDS)
			quot->val[i - j] = q;
		nacc[i--] = 0;
	}

	if (shift) {
		for (k = 0; k < j-1; k++) {
			rem->val[k] = nacc[k] >> shift |
					nacc[k+1] << comp_shift;
		}
		rem->val[k] = nacc[k] >> shift;
	}
	else {
		memcpy(rem->val, nacc, j * sizeof(uint32));
	}

	j = MIN(num->nwords - denom->nwords + 1, MAX_MP_WORDS);
	quot->nwords = num_nonzero_words(quot->val, (uint32)j);
	rem->nwords = num_nonzero_words(rem->val, denom->nwords);
}

/*---------------------------------------------------------------*/
void mp_divrem(mp_t *num, mp_t *denom, mp_t *quot, mp_t *rem) {

	mp_t tmp_quot, tmp_rem;

	if (quot == NULL)
		quot = &tmp_quot;
	if (rem == NULL)
		rem = &tmp_rem;

	if (mp_cmp(num, denom) < 0) {
		/* no division necessary */
		mp_copy(num, rem); 
		mp_clear(quot);
		return;
	}

	mp_clear(rem);	
	if (denom->nwords <= 1) { 		
		/* 1-word division is special-cased */
		rem->val[0] = mp_divrem_1(num, denom->val[0], quot);
		if (rem->val[0] > 0)
			rem->nwords = 1;
		return;
	}

	/* perform the full long division routine */

	mp_clear(quot);
	mp_divrem_core((big_mp_t *)num, denom, quot, rem);

#if 0
	/* an extremely paranoid check for an extremely
	   complex routine */

	if (num->nwords < MAX_MP_WORDS) {
		mp_t test;
		mp_mul(quot, denom, &test);
		mp_add(&test, rem, &test);
		if (mp_cmp(&test, num)) {
			char buf[1000];
			printf("division failed\n");
			printf("%s\n", mp_sprintf(num, 16, buf));
			printf("%s\n", mp_sprintf(denom, 16, buf));
			printf("%s\n", mp_sprintf(quot, 16, buf));
			printf("%s\n", mp_sprintf(rem, 16, buf));
			exit(-1);
		}
	}
#endif
}

/*---------------------------------------------------------------*/
uint32 mp_divrem_1(mp_t *num, uint32 denom, mp_t *quot) {

	int32 i;
	uint32 rem = 0;

	for (i = (int32)num->nwords; i < MAX_MP_WORDS; i++)
		quot->val[i] = 0;
	
	i = num->nwords - 1;
	if (num->val[i] < denom) {
		rem = num->val[i];
		quot->val[i--] = 0;
	}

	while (i >= 0) {

#if (defined(GCC_ASM32X) || defined(GCC_ASM64X)) && !defined(_ICL_WIN_)
		uint32 quot1 = num->val[i];
		ASM_G("divl %2"
			: "+a"(quot1),"+d"(rem)
			: "rm"(denom) : "cc");

		quot->val[i] = quot1;
#else
		uint64 acc = (uint64)rem << 32 | (uint64)num->val[i];
		quot->val[i] = (uint32)(acc / denom);
		rem = (uint32)(acc % denom);
#endif
		i--;
	}

	i = num->nwords;
	quot->nwords = i;
	if (i && quot->val[i-1] == 0)
		quot->nwords--;

	return rem;
}
	
/*---------------------------------------------------------------*/
void mp_modmul(mp_t *a, mp_t *b, mp_t *n, mp_t *res) {

	uint32 i;
	mp_t *small = a;
	mp_t *large = b;
	big_mp_t prod;

	if (small->nwords > large->nwords) {
		small = b;
		large = a;
	}
	if (mp_is_zero(small)) {
		mp_clear(res);
		return;
	}

	memset(&prod, 0, sizeof(big_mp_t));
	for (i = 0; i < small->nwords; i++) 
		mp_addmul_1(large, small->val[i], prod.val + i);
	
	prod.nwords = num_nonzero_words(prod.val, a->nwords + b->nwords);

	mp_mod((mp_t *)&prod, n, res);
}

/*---------------------------------------------------------------*/
uint32 mp_iroot(mp_t *a, uint32 root, mp_t *res) {

	double fp_root;

	if (mp_is_zero(a)) {
		mp_clear(res);
		return 0;
	}

	fp_root = mp_log(a) / M_LN2 / root;

	if (fp_root > 50.0) {
		uint32 num_words, bias;

		num_words = (uint32)((fp_root - 50.0) / 32.0);
		bias = 1 << ((uint32)fp_root - 50 - 32 * num_words);
		fp_root = pow(2.0, fp_root - 32.0 * num_words) + 
						32.0 * bias;
		mp_d2mp(&fp_root, res);
		if (num_words > 0) {
			uint32 i;
			for (i = res->nwords - 1; (int32)i >= 0; i--)
				res->val[num_words + i] = res->val[i];
			res->nwords += num_words;
		}
	}
	else {
		fp_root = pow(2.0, fp_root) + 1.0;
		mp_d2mp(&fp_root, res);
	}

	if (root == 2) {
		while (1) {
			mp_t q;

			mp_div(a, res, &q);
			if (mp_cmp(&q, res) >= 0) {
				mp_mul(res, res, &q);
				return mp_cmp(&q, a);
			}
			mp_add(res, &q, res);
			mp_rshift(res, 1, res);
		}
	}
	else {
		mp_t mp_root;

		mp_clear(&mp_root);
		mp_root.nwords = 1;
		mp_root.val[0] = root - 1;

		while (1) {
			mp_t pow, q;

			mp_pow(res, &mp_root, &pow);
			mp_div(a, &pow, &q);
			if (mp_cmp(&q, res) >= 0) {
				mp_mul(&pow, res, &q);
				return mp_cmp(&q, a);
			}
			mp_mul_1(res, root - 1, res);
			mp_add(res, &q, res);
			mp_divrem_1(res, root, res);
		}
	}
}

/*---------------------------------------------------------------*/
void mp_gcd(mp_t *x, mp_t *y, mp_t *out) {

	mp_t x0, y0;
	mp_t *xptr, *yptr;
	mp_t *tmp;
	mp_t rem;
	int32 sign;

	mp_copy(x, &x0);
	mp_copy(y, &y0);
	if (mp_cmp(x, y) > 0) {
		xptr = &x0;
		yptr = &y0;
	}
	else {
		xptr = &y0;
		yptr = &x0;
	}

	sign = mp_cmp(xptr, yptr);
	mp_clear(out);
	out->nwords = 1;

	while (!mp_is_zero(yptr)) {
		if (xptr->nwords == 1) {
			out->val[0] = mp_gcd_1(xptr->val[0], 
					mp_mod_1(yptr, xptr->val[0]));
			return;
		}
		else if (yptr->nwords == 1) {
			out->val[0] = mp_gcd_1(yptr->val[0],
					mp_mod_1(xptr, yptr->val[0]));
			return;
		}

		if (sign > 0) {
			mp_mod(xptr, yptr, &rem);
			tmp = xptr;
			xptr = yptr;
			yptr = tmp;
		}
		else {
			mp_mod(yptr, xptr, &rem);
		}
		mp_copy(&rem, yptr);
		sign = mp_cmp(xptr, yptr);
	}

	mp_copy(xptr, out);
}

/*---------------------------------------------------------------*/
char * mp_print(mp_t *a, uint32 base, FILE *f, char *scratch) {

	mp_t tmp;
	char *bufptr;
	uint32 next_letter;

	mp_copy(a, &tmp);
	bufptr = scratch + 32 * MAX_MP_WORDS;
	*bufptr = 0;

	do {
		next_letter = mp_divrem_1(&tmp, base, &tmp);
		if (next_letter < 10)
			*(--bufptr) = '0' + next_letter;
		else
			*(--bufptr) = 'a' + (next_letter - 10);
	} while ( !mp_is_zero(&tmp) );
	
	if (f)
		fprintf(f, "%s", bufptr);

	return bufptr;
}

/*---------------------------------------------------------------*/
void mp_str2mp(char *str, mp_t *a, uint32 base) {

	char *str_start, *str_end;
	int32 digit;
	mp_t mult;

	mp_clear(a);

	if (base > 36) {
		return;
	}
	else if (base == 0) {
		if (str[0] == '0' && tolower(str[1]) == 'x') {
			base = 16; 
			str += 2;
		}
		else if (str[0] == '0') {
			base = 8; 
			str++;
		}
		else
			base = 10;
	}

	str_start = str_end = str;
	while (*str_end) {
		digit = tolower(*str_end);
		if (isdigit(digit))
			digit -= '0';
		else if (isalpha(digit))
			digit = digit - 'a' + 10;
		else
			break;

		if (digit >= (int32)base)
			break;
		str_end++;
	}

	mp_clear(&mult);
	mult.nwords = 1;
	mult.val[0] = 1;
	str_end--;

	while( str_end >= str_start && *str_end) {
		digit = tolower(*str_end);
		if (isdigit(digit))
			digit -= '0';
		else if (isalpha(digit))
			digit = digit - 'a' + 10;

		mp_addmul_1(&mult, (uint32)digit, a->val);
		mp_mul_1(&mult, base, &mult);
		str_end--;
	}

	a->nwords = num_nonzero_words(a->val, MAX_MP_WORDS);
}

/*---------------------------------------------------------------*/
uint32 mp_modinv(mp_t *x, mp_t *p, mp_t *res) {

	/* variable names match those of algorithm 2.1.4
	   from Crandall & Pomerance */

	signed_mp_t tmp_a, *a = &tmp_a;
	signed_mp_t tmp_b, *b = &tmp_b;
	signed_mp_t tmp_g, *g = &tmp_g;
	signed_mp_t tmp_u, *u = &tmp_u;
	signed_mp_t tmp_v, *v = &tmp_v;
	signed_mp_t tmp_w, *w = &tmp_w;

	signed_mp_clear(a); a->num.nwords = a->num.val[0] = 1;
	signed_mp_clear(b);
	mp_copy(p, &g->num); g->sign = POSITIVE;
	signed_mp_clear(u);
	signed_mp_clear(v); v->num.nwords = v->num.val[0] = 1;
	mp_copy(x, &w->num); w->sign = POSITIVE;

	while (!mp_is_zero(&w->num)) {

		signed_mp_t q, prod;
		signed_mp_t *t;

		mp_div(&g->num, &w->num, &q.num);
		q.sign = g->sign ^ w->sign;

		signed_mp_mul(&q, u, &prod);
		signed_mp_sub(a, &prod, a);
		t = u; u = a; a = t;

		signed_mp_mul(&q, v, &prod);
		signed_mp_sub(b, &prod, b);
		t = v; v = b; b = t;

		signed_mp_mul(&q, w, &prod);
		signed_mp_sub(g, &prod, g);
		t = w; w = g; g = t;
	}

	if (g->sign != POSITIVE || !mp_is_one(&g->num))
		return 1;

	if (b->sign == NEGATIVE)
		mp_sub(p, &b->num, res);
	else
		mp_copy(&b->num, res);
	return 0;
}

/*---------------------------------------------------------------*/
int32 mp_legendre_1(uint32 a, uint32 p) {

	uint32 x, y, tmp;
	int32 out = 1;

	x = a;
	y = p;
	while (x) {
		while ((x & 1) == 0) {
			x = x / 2;
			if ( (y & 7) == 3 || (y & 7) == 5 )
				out = -out;
		}

		tmp = x;
		x = y;
		y = tmp;

		if ( (x & 3) == 3 && (y & 3) == 3 )
			out = -out;

		x = x % y;
	}
	if (y == 1)
		return out;
	return 0;
}

/*---------------------------------------------------------------*/
int32 mp_legendre(mp_t *a, mp_t *p) {

	mp_t *x, *y, *tmp;
	mp_t tmp_a, tmp_p;
	mp_t rem;
	int32 out = 1;

	mp_copy(a, &tmp_a);
	mp_copy(p, &tmp_p);
	x = &tmp_a;
	y = &tmp_p;

	while (!mp_is_zero(x)) {
		if (x->nwords == 1 && y->nwords == 1)
			return out * mp_legendre_1(x->val[0], y->val[0]);

		while ((x->val[0] & 1) == 0) {
			mp_rshift(x, 1, x);
			if ( (y->val[0] & 7) == 3 || 
			     (y->val[0] & 7) == 5 )
				out = -out;
		}

		tmp = x;
		x = y;
		y = tmp;

		if ( (x->val[0] & 3) == 3 && (y->val[0] & 3) == 3 )
			out = -out;

		mp_mod(x, y, &rem);
		mp_copy(&rem, x);
	}
	if (mp_is_one(y))
		return out;
	return 0;
}

/*---------------------------------------------------------------*/
void mp_expo(mp_t *a, mp_t *b, mp_t *n, mp_t *res) {

	uint32 i, mask;

	mp_clear(res);
	res->nwords = res->val[0] = 1;
	if (mp_is_zero(b))
		return;

	mask = (uint32)(0x80000000);
	i = b->nwords;
	while (mask) {
		if (b->val[i-1] & mask)
			break;
		mask >>= 1;
	}
	
	while (i) {
		mp_modmul(res, res, n, res);

		if (b->val[i-1] & mask)
			mp_modmul(res, a, n, res);

		mask >>= 1;
		if (mask == 0) {
			i--;
			mask = (uint32)(0x80000000);
		}
	}
}
		
/*---------------------------------------------------------------*/
void mp_pow(mp_t *a, mp_t *b, mp_t *res) {

	uint32 i, mask;
	mp_t tmp;

	mp_clear(res);
	res->nwords = res->val[0] = 1;
	if (mp_is_zero(b))
		return;

	mask = (uint32)(0x80000000);
	i = b->nwords;
	while (mask) {
		if (b->val[i-1] & mask)
			break;
		mask >>= 1;
	}
	
	while (i) {
		mp_mul(res, res, &tmp);
		mp_copy(&tmp, res);

		if (b->val[i-1] & mask) {
			mp_mul(res, a, &tmp);
			mp_copy(&tmp, res);
		}

		mask >>= 1;
		if (mask == 0) {
			i--;
			mask = (uint32)(0x80000000);
		}
	}
}
		
/*---------------------------------------------------------------*/
void mp_rand(uint32 bits, mp_t *res, uint32 *seed1, uint32 *seed2) {

	uint32 i;
	uint32 words = (bits + 31) / 32;

	for (i = 0; i < words; i++)
		res->val[i] = get_rand(seed1, seed2);
	for (; i < MAX_MP_WORDS; i++)
		res->val[i] = 0;

	if (bits & 31)
		res->val[words-1] >>= 32 - (bits & 31);
	res->nwords = num_nonzero_words(res->val, words);
}

/*---------------------------------------------------------------*/
uint32 mp_modsqrt_1(uint32 a, uint32 p) {

	uint32 a0 = a;

	if ( (p & 7) == 3 || (p & 7) == 7 ) {
		return mp_expo_1(a0, (p+1)/4, p);
	}
	else if ( (p & 7) == 5 ) {
		uint32 x, y;
		
		if (a0 >= p)
			a0 = a0 % p;
		x = mp_expo_1(a0, (p+3)/8, p);

		if (mp_modmul_1(x, x, p) == a0)
			return x;

		y = mp_expo_1(2, (p-1)/4, p);

		return mp_modmul_1(x, y, p);
	}
	else {
		uint32 d0, d1, a1, s, t, m;
		uint32 i;

		if (a0 == 1)
			return 1;

		for (d0 = 2; d0 < p; d0++) {
			if (mp_legendre_1(d0, p) != -1)
				continue;
	
			t = p - 1;
			s = 0;
			while (!(t & 1)) {
				s++;
				t = t / 2;
			}
	
			a1 = mp_expo_1(a0, t, p);
			d1 = mp_expo_1(d0, t, p);
	
			for (i = 0, m = 0; i < s; i++) {
				uint32 ad;
	
				ad = mp_expo_1(d1, m, p);
				ad = mp_modmul_1(ad, a1, p);
				ad = mp_expo_1(ad, (uint32)(1) << (s-1-i), p);
				if (ad == (p - 1))
					m += (1 << i);
			}
	
			a1 = mp_expo_1(a0, (t+1)/2, p);
			d1 = mp_expo_1(d1, m/2, p);
			return mp_modmul_1(a1, d1, p);
		}
	}

	printf("modsqrt_1 failed\n");
	exit(-1);
}

/*---------------------------------------------------------------*/
void mp_modsqrt(mp_t *a, mp_t *p, mp_t *res,
		uint32 *seed1, uint32 *seed2) {

	mp_t *a0 = a;
	mp_t tmp;

	if ( (p->val[0] & 7) == 3 || (p->val[0] & 7) == 7 ) {
		mp_add_1(p, 1, &tmp);
		mp_rshift(&tmp, 2, &tmp);
		mp_expo(a0, &tmp, p, res);
	}
	else if ( (p->val[0] & 7) == 5 ) {
		mp_t x, base;
		
		mp_add_1(p, 3, &tmp);
		mp_rshift(&tmp, 3, &tmp);
		mp_expo(a0, &tmp, p, res);

		mp_modmul(res, res, p, &x);
		mp_mod(a0, p, &tmp);
		if (!mp_cmp(&x, &tmp))
			return;

		mp_sub_1(p, 1, &tmp);
		mp_rshift(&tmp, 2, &tmp);
		mp_clear(&base);
		base.nwords = 1;
		base.val[0] = 2;
		mp_expo(&base, &tmp, p, &x);
		mp_modmul(&x, res, p, res);
	}
	else {
		mp_t d0, d1, a1, t, m;
		int32 i, j, s;
		uint32 bits = mp_bits(p);

		mp_rand(bits, &d0, seed1, seed2);
		while (mp_legendre(&d0, p) != -1)
			mp_rand(bits, &d0, seed1, seed2);

		mp_sub_1(p, 1, &t);
		s = mp_rjustify(&t, &t);

		mp_expo(a0, &t, p, &a1);
		mp_expo(&d0, &t, p, &d1);
		mp_clear(&m);

		for (i = 0; i < s; i++) {
			mp_t ad;

			mp_expo(&d1, &m, p, &ad);
			mp_modmul(&ad, &a1, p, &ad);
			for (j = 0; j < (s-1-i); j++) {
				mp_modmul(&ad, &ad, p, &ad);
			}
			mp_add_1(&ad, 1, &ad);
			if (!mp_cmp(&ad, p)) {
				m.nwords = i / 32 + 1;
				m.val[i / 32] |= 1 << (i & 31);
			}
		}

		mp_add_1(&t, 1, &t);
		mp_rshift(&t, 1, &t);
		mp_expo(a0, &t, p, &a1);
		mp_rshift(&m, 1, &m);
		mp_expo(&d1, &m, p, &tmp);
		mp_modmul(&a1, &tmp, p, res);
	}
}
		
/*---------------------------------------------------------------*/
void mp_modsqrt2(mp_t *a, mp_t *p, mp_t *res,
		 uint32 *seed1, uint32 *seed2) {

	mp_t p0, p1, r0, r1, r2;

	mp_mod(a, p, &r0);
	mp_modsqrt(a, p, &r1, seed1, seed2);

	mp_mul(&r1, &r1, &r2);
	mp_sub(a, &r2, &r0);
	mp_div(&r0, p, &r2);
	mp_mod(&r2, p, &r0);

	mp_add(&r1, &r1, &r2);
	mp_sub_1(p, 2, &p0);
	mp_expo(&r2, &p0, p, &p1);
	mp_modmul(&p1, &r0, p, &p1);

	mp_mul(p, p, &r0);
	mp_mul(p, &p1, &r2);
	mp_add(&r2, &r1, &r2);
	mp_mod(&r2, &r0, res);

	mp_rshift(&r0, 1, &r1);
	if (mp_cmp(res, &r1) > 0)
		mp_sub(&r0, res, res);
}

/*---------------------------------------------------------------*/
int32 mp_is_prime(mp_t *p, uint32 *seed1, uint32 *seed2) {

	const uint32 factors[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,
				  47,53,59,61,67,71,73,79,83,89,97,101,
				  103,107,109,113,127,131,137,139,149,
				  151,157,163,167,173,179,181,191,193,
				  197,199,211,223,227,229,233,239,241,251};

	uint32 i, j, bits, num_squares;
	mp_t base, tmp, oddpart, p_minus_1;

	for (i = 0; i < sizeof(factors) / sizeof(uint32); i++) {
		if (p->nwords == 1 && p->val[0] == factors[i])
			return 1;

		if (!mp_mod_1(p, factors[i]))
			return 0;
	}

	if (p->nwords == 1 && p->val[0] < 65536)
		return 1;

	mp_sub_1(p, 1, &p_minus_1);
	mp_copy(&p_minus_1, &oddpart);
	bits = mp_bits(p);
	num_squares = mp_rjustify(&oddpart, &oddpart);

	for (i = 0; i < NUM_WITNESSES; i++) {
		mp_rand(bits, &base, seed1, seed2);
		while(mp_is_zero(&base) || mp_is_one(&base) || 
						!mp_cmp(&base, p))
			mp_rand(bits, &base, seed1, seed2);

		if (mp_cmp(&base, p) > 0)
			mp_sub(&base, p, &base);

		mp_expo(&base, &oddpart, p, &tmp);
		if (mp_is_one(&tmp) || !mp_cmp(&tmp, &p_minus_1)) {
		    	continue;
		}

		for (j = 0; j < num_squares - 1; j++) {
			mp_modmul(&tmp, &tmp, p, &tmp);
			if (!mp_cmp(&tmp, &p_minus_1))
				break;
		}

		if (j == num_squares - 1)
			break;
	}

	if (i == NUM_WITNESSES)
		return 1;
	return 0;
}
		
/*---------------------------------------------------------------*/
void mp_random_prime(uint32 bits, mp_t *res,
			uint32 *seed1, uint32 *seed2) {

	mp_rand(bits, res, seed1, seed2);
	res->val[0] |= 1;
	res->nwords = (bits + 31) / 32;
	if (bits & 31)
		res->val[res->nwords - 1] |= 1 << ((bits & 31) - 1);
	else
		res->val[res->nwords - 1] |= 0x80000000;
	
	while (!mp_is_prime(res, seed1, seed2))
		mp_add_1(res, 2, res);
}
		
/*---------------------------------------------------------------*/
uint32 mp_next_prime(mp_t *p, mp_t *res,
		uint32 *seed1, uint32 *seed2) {

	uint32 inc;

	mp_copy(p, res);
	if (res->val[0] & 1) {
		mp_add_1(res, 2, res);
		inc = 2;
	}
	else {
		res->val[0] |= 1;
		inc = 1;
	}

	while (!mp_is_prime(res, seed1, seed2)) {
		mp_add_1(res, 2, res);
		inc += 2;
	}

	return inc;
}
