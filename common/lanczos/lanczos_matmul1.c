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

#include "lanczos.h"

	/* code for handling matrix multiplies when the
	   matrix is in packed format */

/*-------------------------------------------------------------------*/

static void mul_one_med_block(packed_block_t *curr_block,
			uint64 *curr_col, uint64 *curr_b) {

	uint16 *entries = curr_block->d.med_entries;

	while (1) {
		uint64 accum;

#if defined(GCC_ASM64X)
		uint64 i = 0;
		uint64 row = entries[0];
		uint64 count = entries[1];
#else
		uint32 i = 0;
		uint32 row = entries[0];
		uint32 count = entries[1];
#endif

		if (count == 0)
			break;

		/* Unlike the sparse blocks, medium-dense blocks
		   have enough entries that they can be stored in
		   row-major order, with many entries in each row.
		   One iteration of the while loop handles an entire
		   row at a time */

		/* curr_col and curr_b are both cached, so we have to
		   minimize the number of memory accesses and calculate
		   pointers as early as possible */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG)

	#define _txor(k)				\
		"movzwl %%ax, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%eax		\n\t"	\
		"pxor (%2,%%eax,8), %%mm0	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%1,2), %%eax \n\t"	\
		"movzwl %%cx, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%ecx		\n\t"	\
		"pxor (%2,%%ecx,8), %%mm0	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%1,2), %%ecx \n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%3,%1,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%1,2), %%ecx	\n\t"
		"pxor %0, %0			\n\t"
		"pxor %%mm0, %%mm0		\n\t"
		"cmpl $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addl $16, %1			\n\t"
		"cmpl %4, %1			\n\t"
		"jne 0b				\n\t"
		"pxor %%mm0, %0			\n\t"
		"1:				\n\t"

		:"=y"(accum), "+r"(i)
		:"r"(curr_col), "r"(entries),
		 "g"(count & (uint32)(~15))
		:"%eax", "%ecx", "%edx", "%mm0", "memory", "cc");

	#undef _txor

#elif defined(GCC_ASM64X)

    #define _txor(k)				\
		"movzwq %%ax, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rax		\n\t"	\
		"xorq (%2,%%rax,8), %%rsi	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%1,2), %%eax \n\t"	\
		"movzwq %%cx, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rcx		\n\t"	\
		"xorq (%2,%%rcx,8), %%rsi	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%1,2), %%ecx \n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%3,%1,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%1,2), %%ecx	\n\t"
		"xorq %0, %0			\n\t"
		"xorq %%rsi, %%rsi		\n\t"
		"cmpq $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addq $16, %1			\n\t"
		"cmpq %4, %1			\n\t"
		"jne 0b				\n\t"
		"xorq %%rsi, %0			\n\t"
		"1:				\n\t"

		:"=&r"(accum), "+r"(i)
		:"r"(curr_col), "r"(entries), 
		 "g"(count & (uint64)(~15))
		:"%rax", "%rcx", "%rdx", "%rsi", "memory", "cc");

	#undef _txor

#elif defined(MSC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG)

	#define _txor(k)				\
	    ASM_M movzx edx, ax				\
	    ASM_M pxor	mm1, [esi+edx*8]		\
	    ASM_M shr 	eax, 16				\
	    ASM_M pxor 	mm0, [esi+eax*8]		\
	    ASM_M mov 	eax, [2*(k+2+4)+ebx+edi*2]	\
	    ASM_M movzx edx, cx				\
	    ASM_M pxor 	mm1, [esi+edx*8]		\
	    ASM_M shr	ecx, 16				\
	    ASM_M pxor 	mm0, [esi+ecx*8]		\
	    ASM_M mov 	ecx, [2*(k+2+6)+ebx+edi*2]
			
	ASM_M
	{	
		push ebx
		mov esi, curr_col
		mov ebx, entries
		mov edi, i
		mov edx, count
		mov eax, [2*(2+0)+ebx+edi*2]
		mov ecx, [2*(2+2)+ebx+edi*2]
		pxor mm1, mm1
		pxor mm0, mm0
		and edx, ~15
		je L1
		align 16
	L0:	push edx
		_txor(0) _txor(4) _txor(8) _txor(12)
		pop edx
		add edi, 16
		cmp edi, edx
		jne L0
		pxor mm1, mm0
	L1:	movq accum, mm1
		mov i, edi
		pop ebx
	}

	#undef _txor

#else
	accum = 0;
	for (i = 0; i < (count & (uint32)(~15)); i += 16) {
		accum ^= curr_col[entries[i+2+0]] ^
		         curr_col[entries[i+2+1]] ^
		         curr_col[entries[i+2+2]] ^
		         curr_col[entries[i+2+3]] ^
		         curr_col[entries[i+2+4]] ^
		         curr_col[entries[i+2+5]] ^
		         curr_col[entries[i+2+6]] ^
		         curr_col[entries[i+2+7]] ^
		         curr_col[entries[i+2+8]] ^
		         curr_col[entries[i+2+9]] ^
		         curr_col[entries[i+2+10]] ^
		         curr_col[entries[i+2+11]] ^
		         curr_col[entries[i+2+12]] ^
		         curr_col[entries[i+2+13]] ^
		         curr_col[entries[i+2+14]] ^
		         curr_col[entries[i+2+15]];
	}

#endif
		for (; i < count; i++)
			accum ^= curr_col[entries[i+2]];
		curr_b[row] ^= accum;
		entries += count + 2;
	}
}

/*-------------------------------------------------------------------*/
static void mul_one_block(packed_block_t *curr_block,
			uint64 *curr_col, uint64 *curr_b) {

	uint32 i = 0; 
	uint32 num_entries = curr_block->num_entries;
	entry_idx_t *entries = curr_block->d.entries;

	/* unroll by 16, i.e. the number of matrix elements
	   in one cache line (usually). For 32-bit x86, we get
	   a huge performance boost by using either SSE or MMX
	   registers; not because they intrinsically are faster,
	   but because using them cuts the number of memory
	   operations in half, allowing the processor to buffer
	   more xor operations. Also replace two 16-bit loads
	   with a single 32-bit load and extra arithmetic to
	   unpack the array indices */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG)

	#define _txor(x)				\
		"movl 4*" #x "(%1,%0,4), %%eax  \n\t"	\
		"movzwl %%ax, %%ecx             \n\t"	\
		"movq (%2,%%ecx,8), %%mm0       \n\t"	\
		"shrl $16, %%eax                \n\t"	\
		"pxor (%3,%%eax,8), %%mm0       \n\t"	\
		"movq %%mm0, (%2,%%ecx,8)       \n\t"

	ASM_G volatile(
		"cmpl $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)

		"addl $16, %0			\n\t"
		"cmpl %4, %0			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"

		:"+r"(i)
		:"r"(entries), "r"(curr_b), "r"(curr_col), 
		 "g"(num_entries & (uint32)(~15))
		:"%eax", "%ecx", "%mm0", "memory", "cc");

#elif defined(MSC_ASM32A) && defined(HAS_MMX)

	#define _txor(x)				\
		ASM_M mov	eax, [4*x+edi+esi*4]	\
		ASM_M movzx ecx, ax			\
		ASM_M movq 	mm0, [ebx+ecx*8]	\
		ASM_M shr 	eax, 16			\
		ASM_M pxor 	mm0, [edx+eax*8]	\
		ASM_M movq 	[ebx+ecx*8], mm0
	
	ASM_M
	{
		push ebx
		mov esi, i
		mov edi, entries
		mov ebx, curr_b
		mov ecx, num_entries
		mov edx, curr_col
		and ecx, ~15
		je L1
		align 16
	L0:	push ecx
		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)
		pop ecx
		add esi,16
		cmp esi,ecx
		jne L0
	L1:	mov i, esi
		pop ebx
	}

#else
	#define _txor(x) curr_b[entries[i+x].row_off] ^= \
				 curr_col[entries[i+x].col_off]

	for (i = 0; i < (num_entries & (uint32)(~15)); i += 16) {
		#ifdef MANUAL_PREFETCH
		PREFETCH(entries + i + 48);
		#endif

		_txor( 0); _txor( 1); _txor( 2); _txor( 3);
		_txor( 4); _txor( 5); _txor( 6); _txor( 7);
		_txor( 8); _txor( 9); _txor(10); _txor(11);
		_txor(12); _txor(13); _txor(14); _txor(15);
	}
#endif
	#undef _txor

	for (; i < num_entries; i++) {
		curr_b[entries[i].row_off] ^= curr_col[entries[i].col_off];
	}
}

/*-------------------------------------------------------------------*/
void mul_packed_core(void *data, int thread_num)
{
	/* we skip the first matrix row, since it is handled 
	   in the dense function below */

	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;

	uint32 start_block_c = task->block_num * p->superblock_size;
	uint32 num_blocks_c = MIN(p->superblock_size, 
				p->num_block_cols - start_block_c);

	packed_block_t *start_block = p->blocks + start_block_c +
					p->num_block_cols;
	uint64 *x = p->x + start_block_c * p->block_size;
	uint32 i, j;

	for (i = task->task_num; i < p->num_block_rows - 1; 
					i += p->num_threads) {

		packed_block_t *curr_block = start_block + 
					i * p->num_block_cols;
		uint64 *curr_x = x;
		uint32 b_off = i * p->block_size + p->first_block_size;
		uint64 *b = p->b + b_off;

		if (start_block_c == 0) {
			memset(b, 0, MIN(p->block_size, p->nrows - b_off) * 
						sizeof(uint64));
		}

		for (j = 0; j < num_blocks_c; j++) {
			mul_one_block(curr_block, curr_x, b);
			curr_block++;
			curr_x += p->block_size;
		}
	}
}

/*-------------------------------------------------------------------*/
void mul_packed_small_core(void *data, int thread_num)
{
	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;
	thread_data_t *t = p->thread_data + task->task_num;

	uint32 last_task = (task->task_num == p->num_threads - 1);
	uint32 num_blocks = p->num_block_cols / p->num_threads;
	uint32 block_off = num_blocks * task->task_num;
	uint32 off = p->block_size * block_off;
	uint32 vsize = num_blocks * p->block_size;
	uint64 *x = p->x + off;
	uint64 *b = t->tmp_b;
	packed_block_t *curr_block = p->blocks + block_off;
	uint32 i;

	memset(b, 0, MAX(p->first_block_size,
			64 * (1 + (p->num_dense_rows + 63) / 64)) *
			sizeof(uint64));

	if (p->num_threads == 1) {
		vsize = p->ncols;
	}
	else if (last_task) {
		num_blocks = p->num_block_cols - block_off;
		vsize = p->ncols - off;
	}

	for (i = 0; i < num_blocks; i++) {
		mul_one_med_block(curr_block, x, b);
		curr_block++;
		x += p->block_size;
	}

	/* multiply the densest few rows by x (in batches of 64 rows) */

	for (i = 0; i < (p->num_dense_rows + 63) / 64; i++)
		mul_64xN_Nx64(p->dense_blocks[i] + off, 
				p->x + off, b + 64 * i, vsize);
}
