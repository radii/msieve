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

#ifndef _GPU_XFACE_H
#define _GPU_XFACE_H

#if defined(HAVE_CUDA)

#include <util.h>
#include <cuda.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_GPU 4

typedef struct {
	char name[32];
	int32 compute_version_major;
	int32 compute_version_minor;
	int32 clock_speed; /* in kHz */
	int32 num_compute_units;
	int32 constant_mem_size;
	int32 shared_mem_size;
	uint32 global_mem_size;
	int32 registers_per_block;
	int32 max_threads_per_block;
	int32 can_overlap;
	int32 warp_size;
	int32 max_thread_dim[3];
	int32 max_grid_size[3];
	int32 has_timeout;
	CUdevice device_handle;
} gpu_info_t;

typedef struct {
	int32 num_gpu;
	gpu_info_t info[MAX_GPU];
} gpu_config_t;

char * cuGetErrorMessage(CUresult result);

void gpu_init(gpu_config_t *config);

#define CUDA_TRY(func) \
	{ 			 				\
		CUresult status = func;				\
		if (status != CUDA_SUCCESS) {			\
			printf("error (line %d): %s\n", __LINE__,\
				cuGetErrorMessage(status));	\
			exit(-1);				\
		}						\
	}

#define CUDA_ALIGN_PARAM(offset, alignment) \
	(offset) = ((offset) + (alignment) - 1) & ~((alignment) - 1)

#ifdef __cplusplus
}
#endif

#endif /* HAVE_CUDA */

#endif /* !_CUDA_XFACE_H_ */
