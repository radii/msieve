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

#define DEFAULT_DUMP_INTERVAL 2000

#ifdef HAVE_MPI
	#define MPI_NODE_0_START if (obj->mpi_la_row_rank + \
				obj->mpi_la_col_rank == 0) {

	#define MPI_NODE_0_END }
#else
	#define MPI_NODE_0_START /* nothing */
	#define MPI_NODE_0_END /* nothing */
#endif

#define BIT(x) ((uint64)(1) << (x))

static const uint64 bitmask[64] = {
	BIT( 0), BIT( 1), BIT( 2), BIT( 3), BIT( 4), BIT( 5), BIT( 6), BIT( 7),
	BIT( 8), BIT( 9), BIT(10), BIT(11), BIT(12), BIT(13), BIT(14), BIT(15),
	BIT(16), BIT(17), BIT(18), BIT(19), BIT(20), BIT(21), BIT(22), BIT(23),
	BIT(24), BIT(25), BIT(26), BIT(27), BIT(28), BIT(29), BIT(30), BIT(31),
	BIT(32), BIT(33), BIT(34), BIT(35), BIT(36), BIT(37), BIT(38), BIT(39),
	BIT(40), BIT(41), BIT(42), BIT(43), BIT(44), BIT(45), BIT(46), BIT(47),
	BIT(48), BIT(49), BIT(50), BIT(51), BIT(52), BIT(53), BIT(54), BIT(55),
	BIT(56), BIT(57), BIT(58), BIT(59), BIT(60), BIT(61), BIT(62), BIT(63),
};

/*-------------------------------------------------------------------*/
static uint32 form_post_lanczos_matrix(msieve_obj *obj, uint32 *nrows, 
				uint32 *dense_rows_out, 
				uint32 ncols, la_col_t *cols,
				uint64 **post_lanczos_matrix) {

	uint32 i, j, k;
	uint32 num_dense_rows = *dense_rows_out;
	uint32 dense_row_words;
	uint32 new_dense_rows;
	uint32 new_dense_row_words;
	uint32 final_dense_row_words;
	uint64 mask;
	uint64 *submatrix;
	mp_t tmp;

	/* if the matrix is going to have cache blocking applied,
	   proceed but do not form a post-Lanczos matrix if one
	   is not desired. We have to do this because the block
	   matrix multiply expects the number of dense rows to be
	   a multiple of 64.

	   Otherwise, don't do anything if the Lanczos iteration 
	   would finish quickly */

	submatrix = NULL;

#ifdef HAVE_MPI
	/* for a 2-D grid of MPI processes, only the top row
	   of processes construct a post-lanczos matrix */

	if (obj->mpi_la_row_rank != 0)
		return 0;
#endif

	if (*nrows >= MIN_NROWS_TO_PACK ||
	    (POST_LANCZOS_ROWS > 0 && *nrows >= MIN_POST_LANCZOS_DIM)) {

		if (POST_LANCZOS_ROWS > 0) {
			logprintf(obj, "saving the first %u matrix rows "
					"for later\n", POST_LANCZOS_ROWS);
			submatrix = (uint64 *)xmalloc(ncols * sizeof(uint64));
		}
	}
	else {
		return 0;
	}

	mask = (uint64)(-1) >> (64 - POST_LANCZOS_ROWS);
	dense_row_words = (num_dense_rows + 31) / 32;
	mp_clear(&tmp);

	/* we will be removing the first POST_LANCZOS_ROWS rows
	   from the matrix entirely, and packing together the
	   next few rows. The matrix may have dense rows already, 
	   or these rows may be partially or completely sparse, 
	   in which case we'll have to pack them manually. After
	   the post-lanczos rows are removed, the number of dense 
	   rows remaining is a multiple of 64 (minimum of 64) */

	new_dense_rows = MAX(num_dense_rows, POST_LANCZOS_ROWS);
	new_dense_rows += 64 - (new_dense_rows - POST_LANCZOS_ROWS) % 64;
	new_dense_row_words = (new_dense_rows + 31) / 32;
	final_dense_row_words = (new_dense_rows - POST_LANCZOS_ROWS) / 32;

	for (i = 0; i < ncols; i++) {
		uint32 curr_weight = cols[i].weight;
		uint32 *curr_row = cols[i].data;

		/* build up a bitfield of the rows that will be
		   stored in packed format. Start with the rows
		   that are already packed */

		for (j = 0; j < dense_row_words; j++)
			tmp.val[j] = curr_row[curr_weight + j];

		/* add in the rows from the sparse part of the matrix.
		   Entries from these rows are either added to the
		   new dense bitfield, or moved to fill the holes
		   created by packing the first few sparse rows. In 
		   the latter case, the row index must be biased to 
		   reflect the removed rows */

		for (; j < new_dense_row_words; j++)
			tmp.val[j] = 0;

		for (j = k = 0; j < curr_weight; j++) {
			uint32 curr_index = curr_row[j];

			if (curr_index < new_dense_rows)
				tmp.val[curr_index / 32] |= 
						bitmask[curr_index % 32];
			else
				curr_row[k++] = curr_index - POST_LANCZOS_ROWS;
		}

		tmp.nwords = new_dense_row_words;
#if POST_LANCZOS_ROWS > 0
		/* remove the first POST_LANCZOS_ROWS bits from
		   the bitfield */
		submatrix[i] = ((uint64)tmp.val[0] |
				(uint64)tmp.val[1] << 32) & mask;
#endif

		/* move the rest of the bitfield and repack the (hopefully
		   shorter) current column in the heap */
		cols[i].weight = k;
		if (k + final_dense_row_words > 0) {
			cols[i].data = (uint32 *)xrealloc(curr_row, (k + 
						final_dense_row_words) * 
						sizeof(uint32));
			mp_rshift(&tmp, POST_LANCZOS_ROWS, &tmp);
			memcpy(cols[i].data + k, tmp.val, 
						final_dense_row_words * 
						sizeof(uint32));
		}
		else {
			free(cols[i].data);
			cols[i].data = NULL;
		}
	}

	*nrows -= POST_LANCZOS_ROWS;
	*dense_rows_out = new_dense_rows - POST_LANCZOS_ROWS;
	*post_lanczos_matrix = submatrix;
	return (submatrix != NULL);
}

/*-------------------------------------------------------------------*/
static void mul_64x64_64x64(uint64 *a, uint64 *b, uint64 *c ) {

	/* c[][] = x[][] * y[][], where all operands are 64 x 64
	   (i.e. contain 64 words of 64 bits each). The result
	   may overwrite a or b. */

	uint64 ai, bj, accum;
	uint64 tmp[64];
	uint32 i, j;

	for (i = 0; i < 64; i++) {
		j = 0;
		accum = 0;
		ai = a[i];

		while (ai) {
			bj = b[j];
			if (ai & 1)
				accum ^= bj;
			ai >>= 1;
			j++;
		}

		tmp[i] = accum;
	}
	memcpy(c, tmp, sizeof(tmp));
}

/*-----------------------------------------------------------------------*/
static void transpose_64x64(uint64 *a, uint64 *b) {

	uint32 i, j;
	uint64 tmp[64] = {0};

	for (i = 0; i < 64; i++) {
		uint64 word = a[i];
		uint64 mask = bitmask[i];
		for (j = 0; j < 64; j++) {
			if (word & bitmask[j])
				tmp[j] |= mask;
		}
	}
	memcpy(b, tmp, sizeof(tmp));
}

/*-------------------------------------------------------------------*/
static uint32 find_nonsingular_sub(msieve_obj *obj,
				uint64 *t, uint32 *s, 
				uint32 *last_s, uint32 last_dim, 
				uint64 *w) {

	/* given a 64x64 matrix t[][] (i.e. sixty-four
	   64-bit words) and a list of 'last_dim' column 
	   indices enumerated in last_s[]: 
	   
	     - find a submatrix of t that is invertible 
	     - invert it and copy to w[][]
	     - enumerate in s[] the columns represented in w[][] */

	uint32 i, j;
	uint32 dim;
	uint32 cols[64];
	uint64 M[64][2];
	uint64 mask, *row_i, *row_j;
	uint64 m0, m1;

	/* M = [t | I] for I the 64x64 identity matrix */

	for (i = 0; i < 64; i++) {
		M[i][0] = t[i]; 
		M[i][1] = bitmask[i];
	}

	/* put the column indices from last_s[] into the
	   back of cols[], and copy to the beginning of cols[]
	   any column indices not in last_s[] */

	mask = 0;
	for (i = 0; i < last_dim; i++) {
		cols[63 - i] = last_s[i];
		mask |= bitmask[last_s[i]];
	}
	for (i = j = 0; i < 64; i++) {
		if (!(mask & bitmask[i]))
			cols[j++] = i;
	}

	/* compute the inverse of t[][] */

	for (i = dim = 0; i < 64; i++) {
	
		/* find the next pivot row and put in row i */

		mask = bitmask[cols[i]];
		row_i = M[cols[i]];

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		/* if a pivot row was found, eliminate the pivot
		   column from all other rows */

		if (j < 64) {
			for (j = 0; j < 64; j++) {
				row_j = M[cols[j]];
				if ((row_i != row_j) && (row_j[0] & mask)) {
					row_j[0] ^= row_i[0];
					row_j[1] ^= row_i[1];
				}
			}

			/* add the pivot column to the list of 
			   accepted columns */

			s[dim++] = cols[i];
			continue;
		}

		/* otherwise, use the right-hand half of M[]
		   to compensate for the absence of a pivot column */

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		if (j == 64) {
			logprintf(obj, "lanczos error: submatrix "
					"is not invertible\n");
			return 0;
		}
			
		/* eliminate the pivot column from the other rows
		   of the inverse */

		for (j = 0; j < 64; j++) {
			row_j = M[cols[j]];
			if ((row_i != row_j) && (row_j[1] & mask)) {
				row_j[0] ^= row_i[0];
				row_j[1] ^= row_i[1];
			}
		}

		/* wipe out the pivot row */

		row_i[0] = row_i[1] = 0;
	}

	/* the right-hand half of M[] is the desired inverse */
	
	for (i = 0; i < 64; i++) 
		w[i] = M[i][1];

	return dim;
}

/*-----------------------------------------------------------------------*/
static void transpose_vector(uint32 ncols, uint64 *v, uint64 **trans) {

	/* Hideously inefficent routine to transpose a
	   vector v[] of 64-bit words into a 2-D array
	   trans[][] of 64-bit words */

	uint32 i, j;
	uint32 col;
	uint64 mask, word;

	for (i = 0; i < ncols; i++) {
		col = i / 64;
		mask = bitmask[i % 64];
		word = v[i];
		j = 0;
		while (word) {
			if (word & 1)
				trans[j][col] |= mask;
			word = word >> 1;
			j++;
		}
	}
}

/*-----------------------------------------------------------------------*/
static uint32 combine_cols(uint32 ncols, 
			uint64 *x, uint64 *v, 
			uint64 *ax, uint64 *av) {

	/* Once the block Lanczos iteration has finished, 
	   x[] and v[] will contain mostly nullspace vectors
	   between them, as well as possibly some columns
	   that are linear combinations of nullspace vectors.
	   Given vectors ax[] and av[] that are the result of
	   multiplying x[] and v[] by the matrix, this routine 
	   will use Gauss elimination on the columns of [ax | av] 
	   to find all of the linearly dependent columns. The
	   column operations needed to accomplish this are mir-
	   rored in [x | v] and the columns that are independent
	   are skipped. Finally, the dependent columns are copied
	   back into x[] and represent the nullspace vector output
	   of the block Lanczos code. */

	uint32 i, j, k, bitpos, col, col_words;
	uint64 mask;
	uint64 *matrix[128], *amatrix[128], *tmp;

	col_words = (ncols + 63) / 64;

	for (i = 0; i < 128; i++) {
		matrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					     sizeof(uint64));
		amatrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					      sizeof(uint64));
	}

	/* operations on columns can more conveniently become 
	   operations on rows if all the vectors are first
	   transposed */

	transpose_vector(ncols, x, matrix);
	transpose_vector(ncols, ax, amatrix);
	transpose_vector(ncols, v, matrix + 64);
	transpose_vector(ncols, av, amatrix + 64);

	/* Keep eliminating rows until the unprocessed part
	   of amatrix[][] is all zero. The rows where this
	   happens correspond to linearly dependent vectors
	   in the nullspace */

	for (i = bitpos = 0; i < 128 && bitpos < ncols; bitpos++) {

		/* find the next pivot row */

		mask = bitmask[bitpos % 64];
		col = bitpos / 64;
		for (j = i; j < 128; j++) {
			if (amatrix[j][col] & mask) {
				tmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = tmp;
				tmp = amatrix[i];
				amatrix[i] = amatrix[j];
				amatrix[j] = tmp;
				break;
			}
		}
		if (j == 128)
			continue;

		/* a pivot was found; eliminate it from the
		   remaining rows */

		for (j++; j < 128; j++) {
			if (amatrix[j][col] & mask) {

				/* Note that the entire row, *not*
				   just the nonzero part of it, must
				   be eliminated; this is because the
				   corresponding (dense) row of matrix[][]
				   must have the same operation applied */

				for (k = 0; k < col_words; k++) {
					amatrix[j][k] ^= amatrix[i][k];
					matrix[j][k] ^= matrix[i][k];
				}
			}
		}
		i++;
	}

	/* transpose rows i to 64 back into x[]. Pack the
	   dependencies into the low-order bits of x[] */

	for (j = 0; j < ncols; j++) {
		uint64 word = 0;

		col = j / 64;
		mask = bitmask[j % 64];

		for (k = i; k < 64; k++) {
			if (matrix[k][col] & mask)
				word |= bitmask[k - i];
		}
		x[j] = word;
	}

	for (j = 0; j < 128; j++) {
		free(matrix[j]);
		free(amatrix[j]);
	}

	if (i > 64)
		return 0;
	return 64 - i;
}

/*-----------------------------------------------------------------------*/
static void dump_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			uint64 *x, uint64 **vt_v0, uint64 **v, uint64 *v0,
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 max_n, uint32 dim_solved, uint32 iter,
			uint32 s[2][64], uint32 dim1) {

	char buf[256];
	char buf_old[256];
	FILE *dump_fp;
	uint32 status = 1;

#ifdef HAVE_MPI
    
	/* gather x, v[0], v[1], v[2] and v0 into MPI row 0 */

	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : x,
				packed_matrix->nsubcols, 
				MPI_LONG_LONG, x,
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : v[0],
				packed_matrix->nsubcols,  
				MPI_LONG_LONG, v[0],
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : v[1],
				packed_matrix->nsubcols,  
				MPI_LONG_LONG, v[1],
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : v[2],
				packed_matrix->nsubcols,  
				MPI_LONG_LONG, v[2],
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : v0,
				packed_matrix->nsubcols,  
				MPI_LONG_LONG, v0,
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))	
	
	n = packed_matrix->ncols;	
	
	/* pull the full-size vectors into rank 0 */
	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : x,
				n, MPI_LONG_LONG, x,
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : v[0],
				n, MPI_LONG_LONG, v[0],
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : v[1],
				n, MPI_LONG_LONG, v[1],
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : v[2],
				n, MPI_LONG_LONG, v[2],
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : v0, 
				n, MPI_LONG_LONG, v0,
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
	}
#endif

	MPI_NODE_0_START

	sprintf(buf, "%s.chk0", obj->savefile.name);
	sprintf(buf_old, "%s.chk", obj->savefile.name);
	dump_fp = fopen(buf, "wb");
	if (dump_fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	status &= (fwrite(&max_n, sizeof(uint32), (size_t)1, dump_fp)==1);
	status &= (fwrite(&dim_solved, sizeof(uint32), (size_t)1, dump_fp)==1);
	status &= (fwrite(&iter, sizeof(uint32), (size_t)1, dump_fp)==1);

	status &= (fwrite(vt_a_v[1], sizeof(uint64), (size_t)64, dump_fp)==64);
	status &= (fwrite(vt_a2_v[1], sizeof(uint64), (size_t)64, dump_fp)==64);
	status &= (fwrite(winv[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fwrite(winv[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fwrite(vt_v0[0], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fwrite(vt_v0[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fwrite(vt_v0[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fwrite(s[1], sizeof(uint32), (size_t)64, dump_fp) == 64);
	status &= (fwrite(&dim1, sizeof(uint32), (size_t)1, dump_fp) == 1);

	status &= (fwrite(x, sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fwrite(v[0], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fwrite(v[1], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fwrite(v[2], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fwrite(v0, sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	fclose(dump_fp);

	/* only delete an old checkpoint file if the current 
	   checkpoint completed writing. More paranoid: compute a 
	   cryptographic hash of the file and then verify against 
	   the disk image */

	if (status == 0) {
		printf("error: cannot write new checkpoint file\n");
		printf("error: previous checkpoint file not overwritten\n");
		exit(-1);
	}
#if 1
	{ /* let's keep two latest .chk files? */
		char buf_bak[256];
		sprintf(buf_bak, "%s.bak.chk", obj->savefile.name);
		remove(buf_bak);
		rename(buf_old, buf_bak);
	}
#else
	remove(buf_old);
#endif
	if (rename(buf, buf_old)) {
		printf("error: cannot update checkpoint file\n");
		exit(-1);
	}

	MPI_NODE_0_END
}

/*-----------------------------------------------------------------------*/
static void read_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			uint64 *x, uint64 **vt_v0, uint64 **v, uint64 *v0,
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 max_n, uint32 *dim_solved, 
			uint32 *iter, uint32 s[2][64], uint32 *dim1) {

	uint32 read_n;
	uint32 status;
	char buf[256];
	FILE *dump_fp;

	sprintf(buf, "%s.chk", obj->savefile.name);
	dump_fp = fopen(buf, "rb");
	if (dump_fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	status = 1;
	fread(&read_n, sizeof(uint32), (size_t)1, dump_fp);
	if (read_n != max_n) {
		printf("error: unexpected vector size\n");
		exit(-1);
	}
	status &= (fread(dim_solved, sizeof(uint32), (size_t)1, dump_fp) == 1);
	status &= (fread(iter, sizeof(uint32), (size_t)1, dump_fp) == 1);

	status &= (fread(vt_a_v[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_a2_v[1], sizeof(uint64), (size_t)64, dump_fp)==64);
	status &= (fread(winv[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(winv[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[0], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(s[1], sizeof(uint32), (size_t)64, dump_fp) == 64);
	status &= (fread(dim1, sizeof(uint32), (size_t)1, dump_fp) == 1);

	MPI_NODE_0_START
	status &= (fread(x, sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fread(v[0], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fread(v[1], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fread(v[2], sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	status &= (fread(v0, sizeof(uint64), (size_t)max_n, dump_fp)==max_n);
	MPI_NODE_0_END

#ifdef HAVE_MPI
	/* push the full-size vectors to the top grid row */

	if (obj->mpi_ncols > 1 && obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Scatterv(x, packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_col_rank == 0) ?  
					MPI_IN_PLACE : x, 
				n, MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Scatterv(v[0], packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_col_rank == 0) ?  
					MPI_IN_PLACE : v[0],
				n, MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Scatterv(v[1], packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_col_rank == 0) ?  
					MPI_IN_PLACE : v[1], 
				n, MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Scatterv(v[2], packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_col_rank == 0) ?  
					MPI_IN_PLACE : v[2], 
				n, MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Scatterv(v0, packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_col_rank == 0) ?  
					MPI_IN_PLACE : v0, 
				n, MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
	}

	/* duplicate the top grid row across all the grid rows */

	if (obj->mpi_nrows > 1) {
		/* scatter x, v[0], v[1], v[2] and v0 subparts 
		   within each MPI column */

		MPI_TRY(MPI_Scatterv(x, packed_matrix->subcol_counts,
	                         packed_matrix->subcol_offsets, 
	                         MPI_LONG_LONG, 
	                         (obj->mpi_la_row_rank == 0) ?  
	                         MPI_IN_PLACE : x, 
	                         n, MPI_LONG_LONG, 0, 
				 obj->mpi_la_col_grid))
	    
		MPI_TRY(MPI_Scatterv(v[0], packed_matrix->subcol_counts,
	                         packed_matrix->subcol_offsets, 
	                         MPI_LONG_LONG, 
	                         (obj->mpi_la_row_rank == 0) ?  
	                         MPI_IN_PLACE : v[0], 
	                         n, MPI_LONG_LONG, 0, 
				 obj->mpi_la_col_grid))
	    
		MPI_TRY(MPI_Scatterv(v[1], packed_matrix->subcol_counts,
	                         packed_matrix->subcol_offsets, 
	                         MPI_LONG_LONG, 
	                         (obj->mpi_la_row_rank == 0) ?  
	                         MPI_IN_PLACE : v[1], 
	                         n, MPI_LONG_LONG, 0, 
				 obj->mpi_la_col_grid))
	    
		MPI_TRY(MPI_Scatterv(v[2], packed_matrix->subcol_counts,
	                         packed_matrix->subcol_offsets, 
	                         MPI_LONG_LONG, 
	                         (obj->mpi_la_row_rank == 0) ?  
	                         MPI_IN_PLACE : v[2], 
	                         n, MPI_LONG_LONG, 0, 
				 obj->mpi_la_col_grid))
	    
		MPI_TRY(MPI_Scatterv(v0, packed_matrix->subcol_counts,
	                         packed_matrix->subcol_offsets, 
	                         MPI_LONG_LONG, 
	                         (obj->mpi_la_row_rank == 0) ?  
	                         MPI_IN_PLACE : v0, 
	                         n, MPI_LONG_LONG, 0, 
				 obj->mpi_la_col_grid))
	}
#endif

	fclose(dump_fp);
	if (status == 0) {
		printf("error: checkpoint recovery failed\n");
		exit(-1);
	}
}

/*-----------------------------------------------------------------------*/
static void init_lanczos_state(msieve_obj *obj, 
			packed_matrix_t *packed_matrix, uint64 *scratch,
			uint64 *x, uint64 *v0, uint64 **vt_v0, uint64 **v, 
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 s[2][64], uint32 *dim1) {

	uint32 i;

	/* The computed solution 'x' starts off random,
	   and v[0] starts off as B*x. This initial copy
	   of v[0] must be saved off separately */

#ifdef HAVE_MPI
	if (obj->mpi_la_row_rank == 0) {       
#endif
		for (i = 0; i < n; i++) {
			x[i] = v[0][i] = 
			  (uint64)(get_rand(&obj->seed1, &obj->seed2)) << 32 |
		          (uint64)(get_rand(&obj->seed1, &obj->seed2));
		}
#ifdef HAVE_MPI
	}

	/* all nodes work with vectors of size ncols/P */
	n = packed_matrix->nsubcols;
	
	/* scatter x and v[0] subparts within each MPI column */
	MPI_TRY(MPI_Scatterv(x, packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_row_rank == 0) ?  
					MPI_IN_PLACE : x, 
				n, MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Scatterv(v[0], packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets, 
				MPI_LONG_LONG, 
				(obj->mpi_la_row_rank == 0) ?  
					MPI_IN_PLACE : v[0], 
				n, MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
#endif

	mul_sym_NxN_Nx64(packed_matrix, v[0], v[0], scratch);

	memcpy(v0, v[0], n * sizeof(uint64));

	/* Subscripts larger than zero represent past versions of 
	   these quantities, which start off empty (except for the 
	   past version of s[], which contains all the column 
	   indices) */
	   
	memset(v[1], 0, n * sizeof(uint64));
	memset(v[2], 0, n * sizeof(uint64));
	for (i = 0; i < 64; i++) {
		s[1][i] = i;
		vt_a_v[1][i] = 0;
		vt_a2_v[1][i] = 0;
		winv[1][i] = 0;
		winv[2][i] = 0;
		vt_v0[0][i] = 0;
		vt_v0[1][i] = 0;
		vt_v0[2][i] = 0;
	}
	*dim1 = 64;
}

/*-----------------------------------------------------------------------*/
static uint64 * block_lanczos_core(msieve_obj *obj, 
				packed_matrix_t *packed_matrix,
				uint32 *num_deps_found,
				uint64 *post_lanczos_matrix,
				uint32 dump_interval) {
	
	/* Solve Bx = 0 for some nonzero x; the computed
	   solution, containing up to 64 of these nullspace
	   vectors, is returned */

	uint32 n = packed_matrix->ncols;
	uint32 max_n = packed_matrix->max_ncols;
	uint32 alloc_n = MAX(packed_matrix->nrows, packed_matrix->ncols);
	uint64 *vnext, *v[3], *x, *v0;
	uint64 *winv[3], *vt_v0_next;
	uint64 *vt_a_v[2], *vt_a2_v[2], *vt_v0[3];
	uint64 *scratch;
	uint64 *tmp;
	uint32 s[2][64];
	uint64 d[64], e[64], f[64], f2[64];
	uint32 i; 
	uint32 dim0, dim1;
	uint64 mask0, mask1;

	uint32 iter = 0;
	uint32 dim_solved = 0;
	uint32 first_dim_solved = 0;
	uint32 report_interval = 0;
	uint32 check_interval = 0;
	uint32 next_report = 0;
	uint32 log_eta_once = 0;
	uint32 next_check = 0;
	uint32 next_dump = 0;
	time_t first_time;

	if (packed_matrix->num_threads > 1)
		logprintf(obj, "commencing Lanczos iteration (%u threads)\n",
					packed_matrix->num_threads);
	else
		logprintf(obj, "commencing Lanczos iteration\n");

#ifdef HAVE_MPI
	
	/* all nodes work with vectors of size ncols/P */
	n = packed_matrix->nsubcols;
	   
	/* nodes of grid row(col) 0 will gather vectors for node 0 */
	
	if (packed_matrix->mpi_la_row_rank == 0 || 
	    packed_matrix->mpi_la_col_rank == 0)
		alloc_n = MAX(packed_matrix->nrows, packed_matrix->ncols);
	else
		alloc_n = MAX(packed_matrix->nsubrows, packed_matrix->nsubcols);
	
	/* size-n data; the vectors for node 0 are the maximum size
	 but we can limit the vector size for all the other ranks */
	
	MPI_NODE_0_START
	alloc_n = max_n;
	MPI_NODE_0_END
	
	/* we'll need 2 scratch vectors for the matrix multiply */
	scratch = (uint64 *)xmalloc(2 * MAX(packed_matrix->nrows, 
				packed_matrix->ncols) * sizeof(uint64));
#else	
    
	/* size-n data */

	alloc_n = max_n;
	scratch = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
#endif    

	v[0] = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
	v[1] = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
	v[2] = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
	vnext = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
	x = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
	v0 = (uint64 *)xmalloc(alloc_n * sizeof(uint64));
    
	/* 64x64 data */

	winv[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	winv[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	winv[2] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a_v[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a_v[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a2_v[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a2_v[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[2] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0_next = (uint64 *)xmalloc(64 * sizeof(uint64));

	logprintf(obj, "memory use: %.1f MB\n", (double)
			(packed_matrix_sizeof(packed_matrix)) / 1048576);

	/* initialize */

	*num_deps_found = 0;
	iter = 0;
	dim0 = 0;

	if (obj->flags & MSIEVE_FLAG_NFS_LA_RESTART) {
		read_lanczos_state(obj, packed_matrix, 
				x, vt_v0, v, v0, vt_a_v, vt_a2_v,
				winv, packed_matrix->ncols, max_n, 
				&dim_solved, &iter, s, &dim1);
		logprintf(obj, "restarting at iteration %u (dim = %u)\n",
				iter, dim_solved);
	}
	else {
		init_lanczos_state(obj, packed_matrix, scratch, x, 
				v0, vt_v0, v, vt_a_v, vt_a2_v, 
				winv, packed_matrix->ncols, s, &dim1);
	}

	mask1 = 0;
	for (i = 0; i < dim1; i++)
		mask1 |= bitmask[s[1][i]];

	/* determine if the solver will run long enough that
	   it would be worthwhile to report progress */

	first_time = time(NULL);
	if (max_n > 60000 &&
	    obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		  MSIEVE_FLAG_LOG_TO_STDOUT)) {
		if (max_n > 1000000)
			report_interval = 200;
		else if (max_n > 500000)
			report_interval = 500;
		else if (max_n > 100000)
			report_interval = 2000;
		else
			report_interval = 8000;
		first_dim_solved = dim_solved;
		next_report = dim_solved + report_interval;
	}

	if (dump_interval) {
		/* avoid check (at dump) within 4*64 dim + some cushion */
		next_dump = ((dim_solved + 400) / dump_interval + 1) * 
					dump_interval;
		check_interval = 10000;
		/* avoid next_check within 4*64 dim + some cushion */
		next_check = ((dim_solved + 400) / check_interval + 1) * 
					check_interval;
	}

	/* perform the iteration */

	while (1) {
		iter++;

		/* multiply the current v[0] by the matrix and write
		   to vnext */
              
		mul_sym_NxN_Nx64(packed_matrix, v[0], vnext, scratch);
                
		/* compute v0'*A*v0 and (A*v0)'(A*v0) */

		tmul_64xN_Nx64(packed_matrix, v[0], vnext, vt_a_v[0], n);
		tmul_64xN_Nx64(packed_matrix, vnext, vnext, vt_a2_v[0], n);

		/* if the former is orthogonal to itself, then
		   the iteration has finished */

		for (i = 0; i < 64; i++) {
			if (vt_a_v[0][i] != 0)
				break;
		}
		if (i == 64)
			break;

		/* Find the size-'dim0' nonsingular submatrix
		   of v0'*A*v0, invert it, and list the column
		   indices present in the submatrix */

		dim0 = find_nonsingular_sub(obj, vt_a_v[0], s[0], 
					    s[1], dim1, winv[0]);
		if (dim0 == 0)
			break;

		/* mask0 contains one set bit for every column
		   that participates in the inverted submatrix
		   computed above */

		mask0 = 0;
		for (i = 0; i < dim0; i++)
			mask0 |= bitmask[s[0][i]];

		/* The block Lanczos recurrence depends on all columns
		   of v'Av appearing in the current and/or previous iteration. 
		   Verify that condition here
		  
		   Note that the test only applies if this is not the
		   last Lanczos iteration. I'm not sure that this is right, 
		   but the test fails on the last iteration much more often 
		   than would be expected by chance alone, yet ignoring
		   the failure still produces good dependencies 
		   
		   Note that the last iteration typically has dim_solved
		   slightly less than the number of rows, not the number
		   of columns (=n) */
	
		if (dim_solved < packed_matrix->max_nrows - 64) {
			if ((mask0 | mask1) != (uint64)(-1)) {
				logprintf(obj, "lanczos error (dim = %u): "
						"not all columns used\n",
						dim_solved);
				dim0 = 0;
				break;
			}
		}

		/* begin the computation of the next v. First mask
		   off the vectors that are included in this iteration */

		dim_solved += dim0;
		if (mask0 != (uint64)(-1)) {
			for (i = 0; i < n; i++)
				vnext[i] = vnext[i] & mask0;
		}

		/* begin the computation of the next v' * v0. For 
		   the first three iterations, this requires a full 
		   inner product. For all succeeding iterations, the 
		   next v' * v0 is the sum of three 64x64 products 
		   and is stored in vt_v0_next. */

		if (iter < 4) {
			tmul_64xN_Nx64(packed_matrix, v[0], v0, vt_v0[0], n);
		}
		else if (iter == 4) {
			/* v0 is not needed from now on; recycle it 
			   for use as a check vector */
			memset(v0, 0, n * sizeof(uint64));
		}

		/* perform an integrity check on the iteration. This 
		   verifies that the current value of vnext is orthogonal 
		   to the vnext that was computed about check_interval 
		   dimensions ago
		
		   Checks happen on a fixed schedule, as well as 
		   right before a checkpoint file is written */

		if (check_interval && (dim_solved >= next_check ||
#ifndef HAVE_MPI
		    obj->flags & MSIEVE_FLAG_STOP_SIEVING ||
#endif
		    dim_solved >= next_dump)) {

			tmul_64xN_Nx64(packed_matrix, v0, vnext, d, n);
			for (i = 0; i < 64; i++) {
				if (d[i] != (uint64)0) {
					logprintf(obj, "error: corrupt state, "
					       "please restart from "
					       "checkpoint\n");
					printf("\nerror: corrupt state, "
					       "please restart from "
					       "checkpoint\n");
#ifdef HAVE_MPI
					MPI_Abort(MPI_COMM_WORLD, 
							MPI_ERR_ASSERT);
#else
					exit(-1);
#endif
				}
			}
			/* check passed */
			next_check = ((dim_solved + 400) / check_interval + 1) * 
							check_interval;
			memcpy(v0, vnext, n * sizeof(uint64));
		}

		/* compute d, fold it into vnext and update v'*v0 */

		for (i = 0; i < 64; i++)
			d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];

		mul_64x64_64x64(winv[0], d, d);

		for (i = 0; i < 64; i++)
			d[i] = d[i] ^ bitmask[i];

		tmul_Nx64_64x64_acc(packed_matrix, v[0], d, vnext, n);

		transpose_64x64(d, d);
		mul_64x64_64x64(d, vt_v0[0], vt_v0_next);

		/* compute e, fold it into vnext and update v'*v0 */

		mul_64x64_64x64(winv[1], vt_a_v[0], e);

		for (i = 0; i < 64; i++)
			e[i] = e[i] & mask0;

		tmul_Nx64_64x64_acc(packed_matrix, v[1], e, vnext, n);

		transpose_64x64(e, e);
		mul_64x64_64x64(e, vt_v0[1], e);
		for (i = 0; i < 64; i++)
			vt_v0_next[i] = vt_v0_next[i] ^ e[i];

		/* compute f, fold it in. Montgomery shows that 
		   this is unnecessary (f would be zero) if the 
		   previous value of v had full rank */

		if (mask1 != (uint64)(-1)) {
			mul_64x64_64x64(vt_a_v[1], winv[1], f);

			for (i = 0; i < 64; i++)
				f[i] = f[i] ^ bitmask[i];

			mul_64x64_64x64(winv[2], f, f);

			for (i = 0; i < 64; i++)
				f2[i] = ((vt_a2_v[1][i] & mask1) ^ 
					  vt_a_v[1][i]) & mask0;

			mul_64x64_64x64(f, f2, f);

			tmul_Nx64_64x64_acc(packed_matrix, v[2], f, vnext, n);

			transpose_64x64(f, f);
			mul_64x64_64x64(f, vt_v0[2], f);
			for (i = 0; i < 64; i++)
				vt_v0_next[i] = vt_v0_next[i] ^ f[i];
		}

		/* update the computed solution 'x' */

		mul_64x64_64x64(winv[0], vt_v0[0], d);
		tmul_Nx64_64x64_acc(packed_matrix, v[0], d, x, n);

		/* rotate all the variables */

		tmp = v[2]; 
		v[2] = v[1]; 
		v[1] = v[0]; 
		v[0] = vnext; 
		vnext = tmp;
		
		tmp = winv[2]; 
		winv[2] = winv[1]; 
		winv[1] = winv[0]; 
		winv[0] = tmp;
		
		tmp = vt_v0[2]; 
		vt_v0[2] = vt_v0[1]; 
		vt_v0[1] = vt_v0[0]; 
		vt_v0[0] = vt_v0_next; 
		vt_v0_next = tmp;
		
		tmp = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = tmp;
		
		tmp = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = tmp;

		memcpy(s[1], s[0], 64 * sizeof(uint32));
		mask1 = mask0;
		dim1 = dim0;

		MPI_NODE_0_START

		/* possibly print a status update */

		if (report_interval) {
			if (dim_solved >= next_report) {
				time_t curr_time = time(NULL);
				double elapsed = curr_time - first_time;
				uint32 eta = elapsed * (max_n - dim_solved) /
						(dim_solved - first_dim_solved);

				fprintf(stderr, "linear algebra completed %u "
					"of %u dimensions (%1.1f%%, ETA "
					"%dh%2dm)    \r",
					dim_solved, max_n, 100.0 * dim_solved / 
					max_n, eta / 3600, (eta % 3600) / 60);

				/* report the ETA to the logfile once 
				   (wait 6 intervals for a better ETA) */

				if (++log_eta_once == 6) {
					logprintf(obj, "linear algebra at "
						   "%1.1f%%, ETA %dh%2dm\n",
						100.0 * dim_solved / max_n,
						eta / 3600, 
						(eta % 3600) / 60);
				}
				next_report = dim_solved + report_interval;
				fflush(stderr);
			}
		}

		MPI_NODE_0_END

		/* possibly dump a checkpoint file, check for interrupt.

		   Note that MPI cannot reliably dump a checkpoint when
		   interrupted, because multiple MPI processes must 
		   participate in the dump process but there is no real
		   way to signal them all to enter the following at the
		   same time, short of a logical-or of all the obj->flags
		   fields on every iteration. That costs about 3% of the
		   total runtime, and grows with increasing grid size */

		if (dump_interval) {
			if (dim_solved >= next_dump &&
			    dump_interval == DEFAULT_DUMP_INTERVAL) {

				/* the dump interval is the initial one,
				   chosen to accumulate some timing information.
				   Now compute the real dump interval, 
				   calibrated to happen about once per hour.

				   For MPI, the root node computes the
				   dump interval used by everyone */

				MPI_NODE_0_START
				time_t curr_time = time(NULL);
				double elapsed = curr_time - first_time;

				dump_interval = (3600.0 / elapsed) *
					       (dim_solved - first_dim_solved); 
				dump_interval = MAX(dump_interval,
						   DEFAULT_DUMP_INTERVAL + 1);

				/* make the dump interval a multiple of
				   the check interval. If this is not done,
				   eventually we will perform a check and
				   then less than three iterations later
				   will get a dump, which performs another
				   check. The Lanczos recurrence only
				   guarantees that check vectors more than
				   three iterations back will be orthogonal
				   to the current x, so this will cause 
				   spurious failures */

				dump_interval += check_interval -
						dump_interval %
						check_interval;

				logprintf(obj, "checkpointing every %u "
					   "dimensions\n", dump_interval);
				MPI_NODE_0_END
#ifdef HAVE_MPI
				MPI_TRY(MPI_Bcast(&dump_interval, 1, 
						MPI_INT, 0, obj->mpi_la_grid))
#endif
				next_dump = ((dim_solved + 400) / dump_interval + 1) * 
							dump_interval;
				continue;
			}

			if (
#ifndef HAVE_MPI
			    obj->flags & MSIEVE_FLAG_STOP_SIEVING ||
#endif
			    dim_solved >= next_dump) {

				dump_lanczos_state(obj, packed_matrix, 
						   x, vt_v0, v, v0, 
						   vt_a_v, vt_a2_v, winv, 
						   n, max_n, dim_solved, 
						   iter, s, dim1);
				next_dump = ((dim_solved + 400) / dump_interval + 1) * 
							dump_interval;
			}
			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}

	}

	MPI_NODE_0_START
	if (report_interval)
		fprintf(stderr, "\n");
	MPI_NODE_0_END

	logprintf(obj, "lanczos halted after %u iterations (dim = %u)\n", 
					iter, dim_solved);

	/* free unneeded storage */

	free(vnext);
	free(v0);
	free(vt_a_v[0]);
	free(vt_a_v[1]);
	free(vt_a2_v[0]);
	free(vt_a2_v[1]);
	free(winv[0]);
	free(winv[1]);
	free(winv[2]);
	free(vt_v0_next);
	free(vt_v0[0]);
	free(vt_v0[1]);
	free(vt_v0[2]);

	/* if a recoverable failure occurred, start everything
	   over again */

	MPI_NODE_0_START

	if (dim0 == 0 || (obj->flags & MSIEVE_FLAG_STOP_SIEVING)) {
		free(x);
		free(scratch);
		free(v[0]);
		free(v[1]);
		free(v[2]);
		if (dim0 == 0)
			logprintf(obj, "linear algebra failed; retrying...\n");
#ifdef HAVE_MPI
		/* MPI cannot restart gracefully, or shut down
		   gracefully from an interrupt */
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);
#endif
		return NULL;
	}

	MPI_NODE_0_END


	/* convert the output of the iteration to an actual
	   collection of nullspace vectors. Begin by multiplying
	   the output from the iteration by B */

	mul_MxN_Nx64(packed_matrix, x, v[1], scratch);
	mul_MxN_Nx64(packed_matrix, v[0], v[2], scratch);

#ifdef HAVE_MPI
	/* pull the result vectors into rank 0 */
    
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : x,
				packed_matrix->nsubcols, 
				MPI_LONG_LONG, x,
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
				MPI_IN_PLACE : v[0],
				packed_matrix->nsubcols,  
				MPI_LONG_LONG, v[0],
				packed_matrix->subcol_counts,
				packed_matrix->subcol_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_col_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
				MPI_IN_PLACE : v[1],
				packed_matrix->nsubrows, 
				MPI_LONG_LONG, v[1],
				packed_matrix->subrow_counts,
				packed_matrix->subrow_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_row_grid))
	
	MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
				MPI_IN_PLACE : v[2],
				packed_matrix->nsubrows, 
				MPI_LONG_LONG, v[2],
				packed_matrix->subrow_counts,
				packed_matrix->subrow_offsets,
				MPI_LONG_LONG, 0, 
				obj->mpi_la_row_grid))

	n = packed_matrix->ncols;
	
	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : x,
				n, MPI_LONG_LONG, x,
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ? 
						MPI_IN_PLACE : v[0],
				n, MPI_LONG_LONG, v[0],
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))
	}
	if (obj->mpi_la_col_rank == 0) {
		MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
						MPI_IN_PLACE : v[1],
				packed_matrix->nrows, 
				MPI_LONG_LONG, v[1],
				packed_matrix->row_counts,
				packed_matrix->row_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_col_grid))
		MPI_TRY(MPI_Gatherv((obj->mpi_la_row_rank == 0) ? 
						MPI_IN_PLACE : v[2],
				packed_matrix->nrows, 
				MPI_LONG_LONG, v[2],
				packed_matrix->row_counts,
				packed_matrix->row_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_col_grid))
	}
#endif

	MPI_NODE_0_START
        
	/* make sure the last few words of the above matrix products
	   are zero, since the postprocessing will be using them */

	for (i = packed_matrix->max_nrows; 
			i < packed_matrix->max_ncols; i++) {
		v[1][i] = v[2][i] = 0;
	}

	/* if necessary, add in the contribution of the
	   first few rows that were originally in B. We 
	   expect there to be about 64 - POST_LANCZOS_ROWS 
	   bit vectors that are in the nullspace of B and
	   post_lanczos_matrix simultaneously */

	if (post_lanczos_matrix) {
		for (i = 0; i < POST_LANCZOS_ROWS; i++) {
			uint64 accum0 = 0;
			uint64 accum1 = 0;
			uint32 j;
			mask0 = bitmask[i];
			for (j = 0; j < max_n; j++) {
				if (post_lanczos_matrix[j] & mask0) {
					accum0 ^= x[j];
					accum1 ^= v[0][j];
				}
			}
			v[1][i] ^= accum0;
			v[2][i] ^= accum1;
		}
	}

	*num_deps_found = combine_cols(max_n, x, v[0], v[1], v[2]);
	MPI_NODE_0_END

	free(scratch);
	free(v[0]);
	free(v[1]);
	free(v[2]);

	if (*num_deps_found == 0)
		logprintf(obj, "lanczos error: only trivial "
				"dependencies found\n");
	else
		logprintf(obj, "recovered %u nontrivial dependencies\n", 
				*num_deps_found);
	return x;
}

/*-----------------------------------------------------------------------*/
uint64 * block_lanczos(msieve_obj *obj, 
			uint32 nrows, uint32 max_nrows, uint32 start_row,
			uint32 num_dense_rows, 
			uint32 ncols, uint32 max_ncols, uint32 start_col,
			la_col_t *B, uint32 *num_deps_found) {
	
	/* External interface to the linear algebra */

	uint64 *post_lanczos_matrix = NULL;
	uint64 *dependencies;
	packed_matrix_t packed_matrix;
	uint32 dump_interval;
	uint32 have_post_lanczos;
#ifdef HAVE_MPI
	uint32 start_sub;
#endif

	if (max_ncols <= max_nrows) {
		logprintf(obj, "matrix needs more columns than rows; "
                 "try adding 2-3%% more relations\n");
		exit(-1);
	}

	/* optionally remove the densest rows of the matrix, and
	   optionally pack a few more rows into dense format */

	have_post_lanczos = form_post_lanczos_matrix(obj, &nrows,
					&num_dense_rows, ncols, B, 
					&post_lanczos_matrix);
	if (num_dense_rows) {
		logprintf(obj, "matrix includes %u packed rows\n", 
					num_dense_rows);
	}

	memset(&packed_matrix, 0, sizeof(packed_matrix_t));

#ifndef HAVE_MPI
	if (have_post_lanczos)
		max_nrows -= POST_LANCZOS_ROWS;

#else
	/* tell all the MPI processes whether a post lanczos matrix
	   was constructed */

	MPI_TRY(MPI_Bcast(&have_post_lanczos, 1, MPI_INT, 0,
			obj->mpi_la_col_grid))

	if (have_post_lanczos) {
		/* adjust the number of rows to reflect the fact
		   that the matrix is now POST_LANCZOS_ROWS smaller.
		   Fortunately, MPI processes below the top row of
		   the grid have their row numbers relative to offset
		   start_row, so we don't have to adjust all the data */
		
		max_nrows -= POST_LANCZOS_ROWS;
		if (obj->mpi_la_row_rank > 0)
			start_row -= POST_LANCZOS_ROWS;
	}

	/* give the bounds for scatter-gather operations 
	   to the MPI root node */

	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gather(&ncols, 1, MPI_INT, 
				packed_matrix.col_counts, 
				1, MPI_INT, 0, obj->mpi_la_row_grid))
		MPI_TRY(MPI_Gather(&start_col, 1, MPI_INT, 
				packed_matrix.col_offsets, 
				1, MPI_INT, 0, obj->mpi_la_row_grid))
	}

	if (obj->mpi_la_col_rank == 0) {
		MPI_TRY(MPI_Gather(&nrows, 1, MPI_INT, 
				packed_matrix.row_counts, 
				1, MPI_INT, 0, obj->mpi_la_col_grid))
		MPI_TRY(MPI_Gather(&start_row, 1, MPI_INT, 
				packed_matrix.row_offsets, 
				1, MPI_INT, 0, obj->mpi_la_col_grid))
	}
    
	/* figure out the bounds of scatter-gather operations
	   down each MPI column and across each MPI row */

	global_chunk_info(ncols, obj->mpi_nrows, obj->mpi_la_row_rank,
			 &packed_matrix.nsubcols, &start_sub);   
	
	MPI_TRY(MPI_Allgather(&packed_matrix.nsubcols, 1, MPI_INT, 
			packed_matrix.subcol_counts, 
			1, MPI_INT, obj->mpi_la_col_grid))
	MPI_TRY(MPI_Allgather(&start_sub, 1, MPI_INT, 
			packed_matrix.subcol_offsets, 
			1, MPI_INT, obj->mpi_la_col_grid))
	
	global_chunk_info(nrows, obj->mpi_ncols, obj->mpi_la_col_rank,
			&packed_matrix.nsubrows, &start_sub);   
	
	MPI_TRY(MPI_Allgather(&packed_matrix.nsubrows, 1, MPI_INT, 
			packed_matrix.subrow_counts, 
			1, MPI_INT, obj->mpi_la_row_grid))
	MPI_TRY(MPI_Allgather(&start_sub, 1, MPI_INT, 
			packed_matrix.subrow_offsets, 
			1, MPI_INT, obj->mpi_la_row_grid))

	/* if using a post-lanczos matrix, gather the matrix elements
	   at the root node since all of them will be necessary at once */

	if (post_lanczos_matrix != NULL && obj->mpi_la_row_rank == 0) {
		if (obj->mpi_la_col_rank == 0) {
			post_lanczos_matrix = xrealloc(post_lanczos_matrix,
						max_ncols * sizeof(uint64));
		}

		MPI_TRY(MPI_Gatherv((obj->mpi_la_col_rank == 0) ?
					MPI_IN_PLACE : post_lanczos_matrix, 
				ncols, MPI_LONG_LONG, 
				post_lanczos_matrix,
				packed_matrix.col_counts,
				packed_matrix.col_offsets,
				MPI_LONG_LONG, 0, obj->mpi_la_row_grid))

		if (obj->mpi_la_col_rank != 0) {
			free(post_lanczos_matrix);
			post_lanczos_matrix = NULL;
		}
	}
#endif
	if (have_post_lanczos)
		count_matrix_nonzero(obj, nrows, num_dense_rows, ncols, B);

	packed_matrix_init(obj, &packed_matrix, B, 
			   nrows, max_nrows, start_row,
			   ncols, max_ncols, start_col, 
			   num_dense_rows,
#ifdef HAVE_MPI
			   NUM_MEDIUM_ROWS / obj->mpi_nrows
#else
			   NUM_MEDIUM_ROWS
#endif
			   );

	/* set up for writing checkpoint files. This only applies
	   to the largest matrices. The initial dump interval is
	   just to establish timing information */

	dump_interval = 0;
	if (max_nrows > 1000000) {
		dump_interval = DEFAULT_DUMP_INTERVAL;
		obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	}

	/* solve the matrix */

	do {
		dependencies = block_lanczos_core(obj, &packed_matrix,
						num_deps_found,
						post_lanczos_matrix,
						dump_interval);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

	} while (dependencies == NULL);

	if (dump_interval)
		obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	/* note that the following frees any auxiliary packed
	   matrix structures, and also frees the column entries from
	   the input matrix (whether packed or not) */

	packed_matrix_free(&packed_matrix);
	free(post_lanczos_matrix);
	return dependencies;
}
