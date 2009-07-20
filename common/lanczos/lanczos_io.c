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

/*--------------------------------------------------------------------*/
void dump_cycles(msieve_obj *obj, la_col_t *cols, uint32 ncols) {

	uint32 i;
	char buf[256];
	FILE *cycle_fp;

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "wb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: can't open cycle file\n");
		exit(-1);
	}

	fwrite(&ncols, sizeof(uint32), (size_t)1, cycle_fp);

	for (i = 0; i < ncols; i++) {
		la_col_t *c = cols + i;
		uint32 num = c->cycle.num_relations;
		
		fwrite(&num, sizeof(uint32), (size_t)1, cycle_fp);
		fwrite(c->cycle.list, sizeof(uint32), (size_t)num, cycle_fp);
	}
	fclose(cycle_fp);
}

/*--------------------------------------------------------------------*/
void dump_matrix(msieve_obj *obj, 
		uint32 nrows, uint32 num_dense_rows,
		uint32 ncols, la_col_t *cols) {

	uint32 i;
	uint32 dense_row_words;
	char buf[256];
	FILE *matrix_fp;

	dump_cycles(obj, cols, ncols);

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "wb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file\n");
		exit(-1);
	}

	fwrite(&nrows, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&num_dense_rows, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&ncols, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (num_dense_rows + 31) / 32;

	for (i = 0; i < ncols; i++) {
		la_col_t *c = cols + i;
		uint32 num = c->weight + dense_row_words;
		
		fwrite(&c->weight, sizeof(uint32), (size_t)1, matrix_fp);
		fwrite(c->data, sizeof(uint32), (size_t)num, matrix_fp);
	}
	fclose(matrix_fp);
}

/*--------------------------------------------------------------------*/
#define MAX_CYCLE_SIZE 500

void read_cycles(msieve_obj *obj, 
		uint32 *num_cycles_out, 
		la_col_t **cycle_list_out, 
		uint32 dependency) {

	uint32 i;
	uint32 num_cycles;
	uint32 curr_cycle;
	uint32 rel_index[MAX_CYCLE_SIZE];
	char buf[256];
	FILE *cycle_fp;
	FILE *dep_fp = NULL;
	la_col_t *cycle_list;
	uint64 mask = 0;

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "rb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: read_cycles can't open cycle file\n");
		exit(-1);
	}

	if (dependency) {
		sprintf(buf, "%s.dep", obj->savefile.name);
		dep_fp = fopen(buf, "rb");
		if (dep_fp == NULL) {
			logprintf(obj, "error: read_cycles can't "
					"open dependency file\n");
			exit(-1);
		}
		mask = (uint64)1 << (dependency - 1);
	}

	/* read the number of cycles to expect */

	fread(&num_cycles, sizeof(uint32), (size_t)1, cycle_fp);
	cycle_list = (la_col_t *)xcalloc((size_t)num_cycles, sizeof(la_col_t));

	/* read the relation numbers for each cycle */

	for (i = curr_cycle = 0; i < num_cycles; i++) {

		la_col_t *c;
		uint32 num_relations;

		if (fread(&num_relations, sizeof(uint32), 
					(size_t)1, cycle_fp) != 1)
			break;

		if (num_relations > MAX_CYCLE_SIZE) {
			printf("error: cycle too large; corrupt file?\n");
			exit(-1);
		}

		if (fread(rel_index, sizeof(uint32), (size_t)num_relations, 
					cycle_fp) != num_relations)
			break;

		/* all the relation numbers for this cycle
		   have been read; save them and start the
		   count for the next cycle. If reading in 
		   relations to produce a particular dependency
		   from the linear algebra phase, skip any
		   cycles that will not appear in the dependency */

		if (dependency) {
			uint64 curr_dep;

			if (fread(&curr_dep, sizeof(uint64), 
						(size_t)1, dep_fp) == 0) {
				printf("dependency file corrupt\n");
				exit(-1);
			}
			if (!(curr_dep & mask))
				continue;
		}

		c = cycle_list + curr_cycle++;
		c->cycle.num_relations = num_relations;
		c->cycle.list = (uint32 *)xmalloc(num_relations * 
						sizeof(uint32));
		memcpy(c->cycle.list, rel_index, 
				num_relations * sizeof(uint32));
	}
	logprintf(obj, "read %u cycles\n", curr_cycle);
	num_cycles = curr_cycle;
	fclose(cycle_fp);
	if (dep_fp) {
		fclose(dep_fp);
	}
	if (num_cycles == 0) {
		free(cycle_list);
		*num_cycles_out = 0;
		*cycle_list_out = NULL;
		return;
	}

	*num_cycles_out = num_cycles;
	*cycle_list_out = (la_col_t *)xrealloc(cycle_list, 
				num_cycles * sizeof(la_col_t));
}
/*--------------------------------------------------------------------*/
void read_matrix(msieve_obj *obj, 
		uint32 *nrows_out, uint32 *num_dense_rows_out,
		uint32 *ncols_out, la_col_t **cols_out) {

	uint32 i;
	uint32 dense_row_words;
	uint32 ncols;
	la_col_t *cols;
	char buf[256];
	FILE *matrix_fp;

	read_cycles(obj, &ncols, &cols, 0);

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "rb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file\n");
		exit(-1);
	}

	fread(nrows_out, sizeof(uint32), (size_t)1, matrix_fp);
	fread(num_dense_rows_out, sizeof(uint32), (size_t)1, matrix_fp);
	fread(ncols_out, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (*num_dense_rows_out + 31) / 32;
	if (*ncols_out != ncols) {
		printf("error: cycle file not in sync with matrix file\n");
		exit(-1);
	}

	for (i = 0; i < ncols; i++) {
		la_col_t *c = cols + i;
		uint32 num;
		
		fread(&num, sizeof(uint32), (size_t)1, matrix_fp);
		c->weight = num;
		c->data = (uint32 *)xmalloc((num + dense_row_words) * 
					sizeof(uint32));
		fread(c->data, sizeof(uint32), 
				(size_t)(num + dense_row_words), matrix_fp);
	}
	fclose(matrix_fp);
	*cols_out = cols;
}

/*--------------------------------------------------------------------*/
void dump_dependencies(msieve_obj *obj, 
			uint64 *deps, uint32 ncols) {

	char buf[256];
	FILE *deps_fp;

	/* we allow up to 64 dependencies, even though the
	   average case will have (64 - POST_LANCZOS_ROWS) */

	sprintf(buf, "%s.dep", obj->savefile.name);
	deps_fp = fopen(buf, "wb");
	if (deps_fp == NULL) {
		logprintf(obj, "error: can't open deps file\n");
		exit(-1);
	}

	fwrite(deps, sizeof(uint64), (size_t)ncols, deps_fp);
	fclose(deps_fp);
}

