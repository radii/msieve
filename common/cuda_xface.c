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

#include <cuda_xface.h>

#ifdef  HAVE_CUDA

/*------------------------------------------------------------------------*/
char *
cuGetErrorMessage(CUresult result) 
{
	switch (result) {
	case CUDA_SUCCESS: return "CUDA_SUCCESS";
	case CUDA_ERROR_INVALID_VALUE: return "CUDA_ERROR_INVALID_VALUE";
	case CUDA_ERROR_OUT_OF_MEMORY: return "CUDA_ERROR_OUT_OF_MEMORY";
	case CUDA_ERROR_NOT_INITIALIZED: return "CUDA_ERROR_NOT_INITIALIZED";
	case CUDA_ERROR_DEINITIALIZED: return "CUDA_ERROR_DEINITIALIZED";
	case CUDA_ERROR_NO_DEVICE: return "CUDA_ERROR_NO_DEVICE";
	case CUDA_ERROR_INVALID_DEVICE: return "CUDA_ERROR_INVALID_DEVICE";
	case CUDA_ERROR_INVALID_IMAGE: return "CUDA_ERROR_INVALID_IMAGE";
	case CUDA_ERROR_INVALID_CONTEXT: return "CUDA_ERROR_INVALID_CONTEXT";
	case CUDA_ERROR_CONTEXT_ALREADY_CURRENT: return "CUDA_ERROR_CONTEXT_ALREADY_CURRENT";
	case CUDA_ERROR_MAP_FAILED: return "CUDA_ERROR_MAP_FAILED";
	case CUDA_ERROR_UNMAP_FAILED: return "CUDA_ERROR_UNMAP_FAILED";
	case CUDA_ERROR_ARRAY_IS_MAPPED: return "CUDA_ERROR_ARRAY_IS_MAPPED";
	case CUDA_ERROR_ALREADY_MAPPED: return "CUDA_ERROR_ALREADY_MAPPED";
	case CUDA_ERROR_NO_BINARY_FOR_GPU: return "CUDA_ERROR_NO_BINARY_FOR_GPU";
	case CUDA_ERROR_ALREADY_ACQUIRED: return "CUDA_ERROR_ALREADY_ACQUIRED";
	case CUDA_ERROR_NOT_MAPPED: return "CUDA_ERROR_NOT_MAPPED";
	case CUDA_ERROR_INVALID_SOURCE: return "CUDA_ERROR_INVALID_SOURCE";
	case CUDA_ERROR_FILE_NOT_FOUND: return "CUDA_ERROR_FILE_NOT_FOUND";
	case CUDA_ERROR_INVALID_HANDLE: return "CUDA_ERROR_INVALID_HANDLE";
	case CUDA_ERROR_NOT_FOUND: return "CUDA_ERROR_NOT_FOUND";
	case CUDA_ERROR_NOT_READY: return "CUDA_ERROR_NOT_READY";
	case CUDA_ERROR_LAUNCH_FAILED: return "CUDA_ERROR_LAUNCH_FAILED";
	case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES: return "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES";
	case CUDA_ERROR_LAUNCH_TIMEOUT: return "CUDA_ERROR_LAUNCH_TIMEOUT";
	case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING: return "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING";
	case CUDA_ERROR_UNKNOWN: return "CUDA_ERROR_UNKNOWN";
	default: return "unexpected error";
	}
}

/*------------------------------------------------------------------------*/
void
gpu_init(gpu_config_t *config)
{
	int32 i, j;

	/* determine the specifics of CUDA GPUs using the 14 different
	   methods in the Nvidia documentation */

	memset(config, 0, sizeof(gpu_config_t));

	CUDA_TRY(cuInit(0))
	CUDA_TRY(cuDeviceGetCount(&config->num_gpu))
	if (config->num_gpu == 0)
		return;

	for (i = 0; i < (int32)config->num_gpu; i++) {
		CUdevice device;
		CUdevprop prop;
		gpu_info_t *info = config->info + i;

		CUDA_TRY(cuDeviceGet(&device, i))

		info->device_handle = device;

		CUDA_TRY(cuDeviceGetName(info->name,
				sizeof(info->name), device))
		CUDA_TRY(cuDeviceComputeCapability(
				&info->compute_version_major,
				&info->compute_version_minor, device))
		CUDA_TRY(cuDeviceGetProperties(&prop, device))

		info->clock_speed = prop.clockRate;
		info->constant_mem_size = prop.totalConstantMemory;
		info->shared_mem_size = prop.sharedMemPerBlock;
		info->registers_per_block = prop.regsPerBlock;
		info->max_threads_per_block = prop.maxThreadsPerBlock;
		info->warp_size = prop.SIMDWidth;
		for (j = 0; j < 3; j++) {
			info->max_thread_dim[j] = prop.maxThreadsDim[j];
			info->max_grid_size[j] = prop.maxGridSize[j];
		}
		
		CUDA_TRY(cuDeviceTotalMem(
			&info->global_mem_size, device))

		CUDA_TRY(cuDeviceGetAttribute(&info->can_overlap,
				CU_DEVICE_ATTRIBUTE_GPU_OVERLAP, device))

		CUDA_TRY(cuDeviceGetAttribute(&info->num_compute_units,
				CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,
				device))

		CUDA_TRY(cuDeviceGetAttribute(&info->has_timeout,
				CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT,
				device))
	}
}

#endif
