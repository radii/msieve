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

/*-------------------------------------------------------------------*/
static void mul_unpacked(packed_matrix_t *matrix,
			  uint64 *x, uint64 *b) 
{
	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	la_col_t *A = matrix->unpacked_cols;
	uint32 i, j;

	memset(b, 0, ncols * sizeof(uint64));
	
	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		uint64 tmp = x[i];

		for (j = 0; j < col->weight; j++) {
			b[row_entries[j]] ^= tmp;
		}
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			uint64 tmp = x[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] & 
						((uint32)1 << (j % 32))) {
					b[j] ^= tmp;
				}
			}
		}
	}
}

/*-------------------------------------------------------------------*/
static void mul_trans_unpacked(packed_matrix_t *matrix,
				uint64 *x, uint64 *b) 
{
	uint32 ncols = matrix->ncols;
	uint32 num_dense_rows = matrix->num_dense_rows;
	la_col_t *A = matrix->unpacked_cols;
	uint32 i, j;

	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		uint32 *row_entries = col->data;
		uint64 accum = 0;

		for (j = 0; j < col->weight; j++) {
			accum ^= x[row_entries[j]];
		}
		b[i] = accum;
	}

	if (num_dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			uint32 *row_entries = col->data + col->weight;
			uint64 accum = b[i];
	
			for (j = 0; j < num_dense_rows; j++) {
				if (row_entries[j / 32] &
						((uint32)1 << (j % 32))) {
					accum ^= x[j];
				}
			}
			b[i] = accum;
		}
	}
}

/*-------------------------------------------------------------------*/
static void mul_packed(packed_matrix_t *matrix, 
			uint64 *x, uint64 *b) 
{
	uint32 i, j;
	task_control_t task = {NULL, NULL, NULL, NULL};

	matrix->x = x;
	matrix->b = b;

	/* start accumulating the dense matrix multiply results;
	   each thread has scratch space for these, so we don't have
	   to wait for the tasks to finish */

	task.run = mul_packed_small_core;

	for (i = 0; i < matrix->num_threads - 1; i++) {
		task.data = matrix->tasks + i;
		threadpool_add_task(matrix->threadpool, &task, 1);
	}
	mul_packed_small_core(matrix->tasks + i, i);

	/* switch to the sparse blocks */

	task.run = mul_packed_core;

	for (i = 0; i < matrix->num_superblock_cols; i++) {

		la_task_t *t = matrix->tasks;
				
		for (j = 0; j < matrix->num_threads; j++)
			t[j].block_num = i;

		for (j = 0; j < matrix->num_threads - 1; j++) {
			task.data = t + j;
			threadpool_add_task(matrix->threadpool, &task, 1);
		}

		mul_packed_core(t + j, j);
		if (j) {
			threadpool_drain(matrix->threadpool, 1);
		}
	}

	/* xor the small vectors from each thread */

	memcpy(b, matrix->thread_data[0].tmp_b, 
			MAX(matrix->first_block_size,
				64 * ((matrix->num_dense_rows + 63) / 64)) *
			sizeof(uint64));

	for (i = 1; i < matrix->num_threads; i++) {
		accum_xor(b, matrix->thread_data[i].tmp_b, 
				MAX(matrix->first_block_size,
				    64 * ((matrix->num_dense_rows + 63) / 64)));
	}

#if defined(GCC_ASM32A) && defined(HAS_MMX)
	ASM_G volatile ("emms");
#elif defined(MSC_ASM32A) && defined(HAS_MMX)
	ASM_M emms
#endif
}

/*-------------------------------------------------------------------*/
static void mul_trans_packed(packed_matrix_t *matrix, 
			uint64 *x, uint64 *b) 
{
	uint32 i, j;
	task_control_t task = {NULL, NULL, NULL, NULL};

	matrix->x = x;
	matrix->b = b;

	task.run = mul_trans_packed_core;

	for (i = 0; i < matrix->num_superblock_rows; i++) {

		la_task_t *t = matrix->tasks;
				
		for (j = 0; j < matrix->num_threads; j++)
			t[j].block_num = i;

		for (j = 0; j < matrix->num_threads - 1; j++) {
			task.data = t + j;
			threadpool_add_task(matrix->threadpool, &task, 1);
		}

		mul_trans_packed_core(t + j, j);
		if (j) {
			threadpool_drain(matrix->threadpool, 1);
		}
	}

	if (matrix->num_dense_rows) {
		/* add in the dense matrix multiply blocks; these don't 
		   use scratch space, but need all of b to accumulate 
		   results so we have to wait until all tasks finish */

		task.run = mul_trans_packed_small_core;

		for (i = 0; i < matrix->num_threads - 1; i++) {
			task.data = matrix->tasks + i;
			threadpool_add_task(matrix->threadpool, &task, 1);
		}

		mul_trans_packed_small_core(matrix->tasks + i, i);
		if (i) {
			threadpool_drain(matrix->threadpool, 1);
		}
	}

#if defined(GCC_ASM32A) && defined(HAS_MMX)
	ASM_G volatile ("emms");
#elif defined(MSC_ASM32A) && defined(HAS_MMX)
	ASM_M emms
#endif
}

/*-------------------------------------------------------------------*/
static int compare_row_off(const void *x, const void *y) {
	entry_idx_t *xx = (entry_idx_t *)x;
	entry_idx_t *yy = (entry_idx_t *)y;

	if (xx->row_off > yy->row_off)
		return 1;
	if (xx->row_off < yy->row_off)
		return -1;

	return (int)xx->col_off - (int)yy->col_off;
}

/*--------------------------------------------------------------------*/
static void pack_med_block(packed_block_t *b)
{
	uint32 j, k, m;
	uint16 *med_entries;
	entry_idx_t *e;

	/* convert the first block in the stripe to a somewhat-
	   compressed format. Entries in this first block are stored 
	   by row, and all rows are concatenated into a single 
	   16-bit array */

	e = b->d.entries;
	qsort(e, (size_t)b->num_entries, 
			sizeof(entry_idx_t), compare_row_off);
	for (j = k = 1; j < b->num_entries; j++) {
		if (e[j].row_off != e[j-1].row_off)
			k++;
	}

	/* we need a 16-bit word for each element and two more
	   16-bit words at the start of each of the k packed
	   arrays making up med_entries. The first extra word
	   gives the row number and the second gives the number
	   of entries in that row. We also need a few extra words 
	   at the array end because the multiply code uses a 
	   software pipeline and would fetch off the end of 
	   med_entries otherwise */

	med_entries = (uint16 *)xmalloc((b->num_entries + 
					2 * k + 8) * sizeof(uint16));
	j = k = 0;
	while (j < b->num_entries) {
		for (m = 0; j + m < b->num_entries; m++) {
			if (m > 0 && e[j+m].row_off != e[j+m-1].row_off)
				break;
			med_entries[k+m+2] = e[j+m].col_off;
		}
		med_entries[k] = e[j].row_off;
		med_entries[k+1] = m;
		j += m;
		k += m + 2;
	}
	med_entries[k] = med_entries[k+1] = 0;
	free(b->d.entries);
	b->d.med_entries = med_entries;
}

/*--------------------------------------------------------------------*/
static void pack_matrix_core(packed_matrix_t *p, la_col_t *A)
{
	uint32 i, j, k;
	uint32 dense_row_blocks;
	packed_block_t *curr_stripe;

	uint32 ncols = p->ncols;
	uint32 block_size = p->block_size;
	uint32 num_block_rows = p->num_block_rows;
	uint32 num_block_cols = p->num_block_cols;
	uint32 num_dense_rows = p->num_dense_rows;
	uint32 first_block_size = p->first_block_size;

	/* pack the dense rows 64 at a time */

	dense_row_blocks = (num_dense_rows + 63) / 64;
	if (dense_row_blocks) {
		p->dense_blocks = (uint64 **)xmalloc(dense_row_blocks *
						sizeof(uint64 *));
		for (i = 0; i < dense_row_blocks; i++) {
			p->dense_blocks[i] = (uint64 *)xmalloc(ncols *
							sizeof(uint64));
		}

		for (i = 0; i < ncols; i++) {
			la_col_t *c = A + i;
			uint32 *src = c->data + c->weight;
			for (j = 0; j < dense_row_blocks; j++) {
				p->dense_blocks[j][i] = 
						(uint64)src[2 * j + 1] << 32 |
						(uint64)src[2 * j];
			}
		}
	}

	/* allocate blocks in row-major order; a 'stripe' is
	   a vertical column of blocks. The first block in each
	   column has first_block_size rows instead of block_size */

	p->blocks = curr_stripe = (packed_block_t *)xcalloc(
						(size_t)num_block_rows *
						        num_block_cols,
						sizeof(packed_block_t));

	/* we convert the sparse part of the matrix to packed
	   format one stripe at a time. This limits the worst-
	   case memory use of the packing process */

	for (i = 0; i < num_block_cols; i++, curr_stripe++) {

		uint32 curr_cols = MIN(block_size, ncols - i * block_size);
		packed_block_t *b;

		/* count the number of nonzero entries in each block */

		for (j = 0; j < curr_cols; j++) {
			la_col_t *c = A + i * block_size + j;
			uint32 limit = first_block_size;

			for (k = 0, b = curr_stripe; k < c->weight; k++) {
				uint32 index = c->data[k];

				while (index >= limit) {
					b += num_block_cols;
					limit += block_size;
				}
				b->num_entries++;
			}
		}

		/* concatenate the nonzero elements of the matrix
		   columns corresponding to this stripe.
		   
		   We technically can combine the previous pass through
		   the columns with this pass, but on some versions of
		   libc the number of reallocations causes an incredible
		   slowdown */

		for (j = 0, b = curr_stripe; j < num_block_rows; 
						j++, b += num_block_cols) {
			b->d.entries = (entry_idx_t *)xmalloc(
						b->num_entries *
						sizeof(entry_idx_t));
			b->num_entries = 0;
		}

		for (j = 0; j < curr_cols; j++) {
			la_col_t *c = A + i * block_size + j;
			uint32 limit = first_block_size;
			uint32 start_row = 0;

			for (k = 0, b = curr_stripe; k < c->weight; k++) {
				entry_idx_t *e;
				uint32 index = c->data[k];

				while (index >= limit) {
					b += num_block_cols;
					start_row = limit;
					limit += block_size;
				}

				e = b->d.entries + b->num_entries++;
				e->row_off = index - start_row;
				e->col_off = j;
			}

			free(c->data);
			c->data = NULL;
		}

		pack_med_block(curr_stripe);
	}
}

/*--------------------------------------------------------------------*/
static void matrix_thread_init(void *data, int thread_num) {

	packed_matrix_t *p = (packed_matrix_t *)data;
	thread_data_t *t = p->thread_data + thread_num;

	/* we use this scratch vector for both matrix multiplies
	   and vector-vector operations; it has to be large enough
	   to support both. Note that first_block_size is split across
	   MPI rows, so it is conceivable with enough MPI processes that
	   the MAX() is necessary */

	t->tmp_b = (uint64 *)xmalloc(MAX(p->first_block_size,
				64 * (1 + (p->num_dense_rows + 63) / 64)) *
				sizeof(uint64));
}

/*-------------------------------------------------------------------*/
static void matrix_thread_free(void *data, int thread_num) {

	packed_matrix_t *p = (packed_matrix_t *)data;
	thread_data_t *t = p->thread_data + thread_num;

	free(t->tmp_b);
}

/*-------------------------------------------------------------------*/
void packed_matrix_init(msieve_obj *obj,
			packed_matrix_t *p, la_col_t *A,
			uint32 nrows, uint32 max_nrows, uint32 start_row, 
			uint32 ncols, uint32 max_ncols, uint32 start_col, 
			uint32 num_dense_rows, uint32 first_block_size) {

	uint32 i, j;
	uint32 block_size;
	uint32 superblock_size;
	uint32 num_threads;
	thread_control_t control;

	/* initialize */

	p->unpacked_cols = A;
	p->nrows = nrows;
	p->max_nrows = max_nrows;
	p->start_row = start_row;
	p->ncols = ncols;
	p->max_ncols = max_ncols;
	p->start_col = start_col;
	p->num_dense_rows = num_dense_rows;
	p->num_threads = 1;
	p->first_block_size = first_block_size; /* needed for thread pool init */
#ifdef HAVE_MPI
	p->mpi_size = obj->mpi_size;
	p->mpi_nrows = obj->mpi_nrows;
	p->mpi_ncols = obj->mpi_ncols;
	p->mpi_la_row_rank = obj->mpi_la_row_rank;
	p->mpi_la_col_rank = obj->mpi_la_col_rank;
	p->mpi_la_row_grid = obj->mpi_la_row_grid;
	p->mpi_la_col_grid = obj->mpi_la_col_grid;
#endif

	/* determine the number of threads to use */

	num_threads = obj->num_threads;
	if (num_threads < 2 || max_nrows < MIN_NROWS_TO_THREAD)
		num_threads = 1;

	p->num_threads = num_threads = MIN(num_threads, MAX_THREADS);

	/* start the thread pool */

	control.init = matrix_thread_init;
	control.shutdown = matrix_thread_free;
	control.data = p;

	if (num_threads > 1) {
		p->threadpool = threadpool_init(num_threads - 1, 
						200, &control);
	}
	matrix_thread_init(p, num_threads - 1);

	/* pre-generate the structures to drive the thread pool;
	   do this even for single-threaded runs */

	p->tasks = (la_task_t *)xmalloc(sizeof(la_task_t) * 
					p->num_threads);

	for (i = 0; i < p->num_threads; i++) {
		p->tasks[i].matrix = p;
		p->tasks[i].task_num = i;
	}

	if (max_nrows <= MIN_NROWS_TO_PACK)
		return;

	/* determine the block sizes. We assume that the largest
	   cache in the system is unified and shared across all
	   threads. When performing matrix multiplies 'A*x=b', 
	   we choose the block size small enough so that one block
	   of b fits in L1 cache, and choose the superblock size
	   to be small enough so that a superblock's worth of x
	   or b takes up 3/4 of the largest cache in the system.
	   
	   Making the block size too small increases memory use
	   and puts more pressure on the larger caches, while
	   making the superblock too small reduces the effectiveness
	   of L1 cache and increases the synchronization overhead
	   in multithreaded runs */

	block_size = 8192;
	superblock_size = 3 * obj->cache_size2 / (4 * sizeof(uint64));

	/* possibly override from the command line */

	if (obj->nfs_args != NULL) {

		const char *tmp;

		tmp = strstr(obj->nfs_args, "la_block=");
		if (tmp != NULL)
			block_size = atoi(tmp + 9);

		tmp = strstr(obj->nfs_args, "la_superblock=");
		if (tmp != NULL)
			superblock_size = atoi(tmp + 14);
	}

	logprintf(obj, "using block size %u and superblock size %u for "
			"processor cache size %u kB\n", 
				block_size, superblock_size,
				obj->cache_size2 / 1024);

	p->unpacked_cols = NULL;

	p->block_size = block_size;
	p->num_block_cols = (ncols + block_size - 1) / block_size;
	p->num_block_rows = 1 + (nrows - first_block_size + 
				block_size - 1) / block_size;

	p->superblock_size = (superblock_size + block_size - 1) / block_size;
	p->num_superblock_cols = (p->num_block_cols + p->superblock_size - 1) / 
					p->superblock_size;
	p->num_superblock_rows = (p->num_block_rows - 1 + p->superblock_size - 1) / 
					p->superblock_size;

	/* do the core work of packing the matrix */

	pack_matrix_core(p, A);
}

/*-------------------------------------------------------------------*/
void packed_matrix_free(packed_matrix_t *p) {

	uint32 i;

	if (p->unpacked_cols) {
		la_col_t *A = p->unpacked_cols;
		for (i = 0; i < p->ncols; i++) {
			free(A[i].data);
			A[i].data = NULL;
		}
	}
	else {
		for (i = 0; i < (p->num_dense_rows + 63) / 64; i++)
			free(p->dense_blocks[i]);

		for (i = 0; i < p->num_block_rows * p->num_block_cols; i++) 
			free(p->blocks[i].d.entries);

		free(p->blocks);
	}

	if (p->num_threads > 1) {
		threadpool_drain(p->threadpool, 1);
		threadpool_free(p->threadpool);
	}
	matrix_thread_free(p, p->num_threads - 1);
	free(p->tasks);
}

/*-------------------------------------------------------------------*/
size_t packed_matrix_sizeof(packed_matrix_t *p) {

	uint32 i, j;
	size_t mem_use;

	/* account for the vectors used in the lanczos iteration */

	if (p->start_row + p->start_col == 0)
		mem_use = 7 * p->max_ncols;
	else
		mem_use = 7 * MAX(p->nrows, p->ncols);

	/* and for the matrix */

	if (p->unpacked_cols) {
		la_col_t *A = p->unpacked_cols;
		mem_use += p->ncols * (sizeof(la_col_t) +
				(p->num_dense_rows + 31) / 32);
		for (i = 0; i < p->ncols; i++) {
			mem_use += A[i].weight * sizeof(uint32);
		}
	}
	else {
		uint32 num_blocks = p->num_block_rows * 
					p->num_block_cols;

		mem_use += sizeof(uint64) * p->num_threads * p->first_block_size;

		mem_use += sizeof(packed_block_t) * num_blocks;

		mem_use += p->ncols * sizeof(uint64) *
				((p->num_dense_rows + 63) / 64);

		for (j = 0; j < num_blocks; j++) {
			packed_block_t *b = p->blocks + j;

			if (j < p->num_block_cols) {
				mem_use += (b->num_entries + 
					    2 * p->first_block_size) * 
						sizeof(uint16);
			}
			else {
				mem_use += b->num_entries *
						sizeof(entry_idx_t);
			}
		}
	}
	return mem_use;
}

/*-------------------------------------------------------------------*/
void mul_MxN_Nx64(packed_matrix_t *A, uint64 *x, 
			uint64 *b, uint64 *scratch) {
    
	/* Multiply the vector x[] by the matrix A (stored
	   columnwise) and put the result in b[]. The MPI 
	   version needs a scratch array because MPI reduction
	   operations apparently cannot be performed in-place */

#ifdef HAVE_MPI
	uint64 *scratch2 = scratch + MAX(A->ncols, A->nrows);

	if (A->mpi_size <= 1) {
#endif
		if (A->unpacked_cols)
			mul_unpacked(A, x, b);
		else
			mul_packed(A, x, b);
#ifdef HAVE_MPI
		return;
	}
    
	/* make each MPI column gather its own part of x */
	
	global_allgather(x, scratch, A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
		
	mul_packed(A, scratch, scratch2);
	
	/* make each MPI row combine and scatter its own part of A^T * A*x */
	
	global_xor_scatter(scratch2, b, scratch, A->nrows, A->mpi_ncols,
			A->mpi_la_col_rank, A->mpi_la_row_grid);

#endif
}

/*-------------------------------------------------------------------*/
void mul_sym_NxN_Nx64(packed_matrix_t *A, uint64 *x, 
			uint64 *b, uint64 *scratch) {

	/* Multiply x by A and write to scratch, then
	   multiply scratch by the transpose of A and
	   write to b. x may alias b, but the two must
	   be distinct from scratch */

#ifdef HAVE_MPI
	uint64 *scratch2 = scratch + MAX(A->ncols, A->nrows);
        
	if (A->mpi_size <= 1) {
#endif
		if (A->unpacked_cols) {
			mul_unpacked(A, x, scratch);
			mul_trans_unpacked(A, scratch, b);
		}
		else {
			mul_packed(A, x, scratch);
			mul_trans_packed(A, scratch, b);
		}
#ifdef HAVE_MPI
		return;
	}
    
	/* make each MPI column gather its own part of x */
	 
	global_allgather(x, scratch, A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
	
	mul_packed(A, scratch, scratch2);
		
	/* make each MPI row combine its own part of A*x */
	
	global_xor(scratch2, scratch, A->nrows, A->mpi_ncols,
			   A->mpi_la_col_rank, A->mpi_la_row_grid);
		
	mul_trans_packed(A, scratch, scratch2);
		
	/* make each MPI row combine and scatter its own part of A^T * A*x */
		
	global_xor_scatter(scratch2, b, scratch,  A->ncols, A->mpi_nrows, 
			A->mpi_la_row_rank, A->mpi_la_col_grid);
#endif
}
