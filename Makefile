# --------------------------------------------------------------------
# This source distribution is placed in the public domain by its author,
# Jason Papadopoulos. You may use it for any purpose, free of charge,
# without having to notify anyone. I disclaim any responsibility for any
# errors.
# 
# Optionally, please be nice and tell me if you find this source to be
# useful. Again optionally, if you add to the functionality present here
# please consider making those additions public too, so that others may 
# benefit from your work.	
#
#  $Id$
# --------------------------------------------------------------------

# override from command line
WIN = 0
WIN64 = 0

# gcc with basic optimization (-march flag could
# get overridden by architecture-specific builds)
CC = gcc
WARN_FLAGS = -Wall -W
OPT_FLAGS = -O3 -fomit-frame-pointer -march=native \
	    -D_FILE_OFFSET_BITS=64 -DNDEBUG -D_LARGEFILE64_SOURCE

# use := instead of = so we only run the following once
SVN_VERSION := $(shell svnversion .)
ifeq ($(SVN_VERSION),)
	SVN_VERSION := unknown
endif

CFLAGS = $(OPT_FLAGS) $(MACHINE_FLAGS) $(WARN_FLAGS) \
	 	-DMSIEVE_SVN_VERSION="\"$(SVN_VERSION)\"" \
		-I. -Iaprcl -Iinclude -Ignfs -Ignfs/poly -Ignfs/poly/stage1

# tweak the compile flags

ifeq ($(ECM),1)
	CFLAGS += -DHAVE_GMP_ECM
	LIBS += -lecm
endif
ifeq ($(WIN),1)

else
	LIBS += -ldl
endif
ifeq ($(CUDA),1)

ifeq ($(WIN),1)
	CUDA_ROOT = $(shell echo $$CUDA_PATH)
	NVCC = "$(CUDA_ROOT)/bin/nvcc"

ifeq ($(WIN64),1)
	CUDA_LIBS = "$(CUDA_ROOT)/lib/x64/cuda.lib"
else
	CUDA_LIBS = "$(CUDA_ROOT)/lib/win32/cuda.lib"
endif

else
	NVCC = "$(shell which nvcc)"
	CUDA_ROOT = $(shell dirname $(NVCC))/../
	CUDA_LIBS = -lcuda
endif
	CFLAGS += -I"$(CUDA_ROOT)/include" -Icub -DHAVE_CUDA
	LIBS += $(CUDA_LIBS)
endif
ifeq ($(MPI),1)
	CC = mpicc
	CFLAGS += -DHAVE_MPI
endif
ifeq ($(BOINC),1)
	# fill in as appropriate
	BOINC_INC_DIR = .
	BOINC_LIB_DIR = .
	CFLAGS += -I$(BOINC_INC_DIR) -DHAVE_BOINC
	LIBS += -L$(BOINC_LIB_DIR) -lboinc_api -lboinc
endif
ifeq ($(NO_ZLIB),1)
	CFLAGS += -DNO_ZLIB
else
	LIBS += -lz
endif


# Note to MinGW users: the library does not use pthreads calls in
# win32 or win64, so it's safe to pull libpthread into the link line.
# Of course this does mean you have to install the minGW pthreads bundle...

LIBS += -lgmp -lm -lpthread

#---------------------------------- Generic file lists -------------------

COMMON_HDR = \
	aprcl/mpz_aprcl32.h \
	common/lanczos/lanczos.h \
	common/filter/filter.h \
	common/filter/filter_priv.h \
	common/filter/merge_util.h \
	include/batch_factor.h \
	include/common.h \
	include/cuda_xface.h \
	include/dd.h \
	include/ddcomplex.h \
	include/gmp_xface.h \
	include/integrate.h \
	include/msieve.h \
	include/mp.h \
	include/polyroot.h \
	include/thread.h \
	include/util.h

COMMON_SRCS = \
	aprcl/mpz_aprcl32.c \
	common/filter/clique.c \
	common/filter/filter.c \
	common/filter/merge.c \
	common/filter/merge_post.c \
	common/filter/merge_pre.c \
	common/filter/merge_util.c \
	common/filter/singleton.c \
	common/lanczos/lanczos.c \
	common/lanczos/lanczos_io.c \
	common/lanczos/lanczos_matmul0.c \
	common/lanczos/lanczos_matmul1.c \
	common/lanczos/lanczos_matmul2.c \
	common/lanczos/lanczos_pre.c \
	common/lanczos/lanczos_vv.c \
	common/lanczos/matmul_util.c \
	common/smallfact/gmp_ecm.c \
	common/smallfact/smallfact.c \
	common/smallfact/squfof.c \
	common/smallfact/tinyqs.c \
	common/batch_factor.c \
	common/cuda_xface.c \
	common/dickman.c \
	common/driver.c \
	common/expr_eval.c \
	common/hashtable.c \
	common/integrate.c \
	common/minimize.c \
	common/minimize_global.c \
	common/mp.c \
	common/polyroot.c \
	common/prime_delta.c \
	common/prime_sieve.c \
	common/savefile.c \
	common/strtoll.c \
	common/thread.c \
	common/util.c

COMMON_OBJS = $(COMMON_SRCS:.c=.o)

#---------------------------------- QS file lists -------------------------

QS_HDR = mpqs/mpqs.h

QS_SRCS = \
	mpqs/gf2.c \
	mpqs/mpqs.c \
	mpqs/poly.c \
	mpqs/relation.c \
	mpqs/sieve.c \
	mpqs/sieve_core.c \
	mpqs/sqrt.c

QS_OBJS = \
	mpqs/gf2.qo \
	mpqs/mpqs.qo \
	mpqs/poly.qo \
	mpqs/relation.qo \
	mpqs/sieve.qo \
	mpqs/sqrt.qo \
	mpqs/sieve_core_generic_32k.qo \
	mpqs/sieve_core_generic_64k.qo

#---------------------------------- GPU file lists -------------------------

GPU_OBJS = \
	stage1_core_sm20.ptx \
	stage1_core_sm30.ptx \
	stage1_core_sm35.ptx \
	cub/built

#---------------------------------- NFS file lists -------------------------

NFS_HDR = \
	gnfs/filter/filter.h \
	gnfs/poly/poly.h \
	gnfs/poly/poly_skew.h \
	gnfs/poly/stage1/stage1.h \
	gnfs/poly/stage2/stage2.h \
	gnfs/sieve/sieve.h \
	gnfs/sqrt/sqrt.h \
	gnfs/gnfs.h

NFS_GPU_HDR = \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core.cu \
	gnfs/poly/stage1/stage1_core_gpu/cuda_intrinsics.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core.h

NFS_NOGPU_HDR = \
	gnfs/poly/stage1/cpu_intrinsics.h

NFS_SRCS = \
	gnfs/poly/poly.c \
	gnfs/poly/poly_param.c \
	gnfs/poly/poly_skew.c \
	gnfs/poly/polyutil.c \
	gnfs/poly/root_score.c \
	gnfs/poly/size_score.c \
	gnfs/poly/stage1/stage1.c \
	gnfs/poly/stage1/stage1_roots.c \
	gnfs/poly/stage2/optimize.c \
	gnfs/poly/stage2/optimize_deg6.c \
	gnfs/poly/stage2/root_sieve.c \
	gnfs/poly/stage2/root_sieve_deg45_x.c \
	gnfs/poly/stage2/root_sieve_deg5_xy.c \
	gnfs/poly/stage2/root_sieve_deg6_x.c \
	gnfs/poly/stage2/root_sieve_deg6_xy.c \
	gnfs/poly/stage2/root_sieve_deg6_xyz.c \
	gnfs/poly/stage2/root_sieve_line.c \
	gnfs/poly/stage2/root_sieve_util.c \
	gnfs/poly/stage2/stage2.c \
	gnfs/filter/duplicate.c \
	gnfs/filter/filter.c \
	gnfs/filter/singleton.c \
	gnfs/sieve/sieve_line.c \
	gnfs/sieve/sieve_util.c \
	gnfs/sqrt/sqrt.c \
	gnfs/sqrt/sqrt_a.c \
	gnfs/fb.c \
	gnfs/ffpoly.c \
	gnfs/gf2.c \
	gnfs/gnfs.c \
	gnfs/relation.c

NFS_OBJS = $(NFS_SRCS:.c=.no)

NFS_GPU_SRCS = \
	gnfs/poly/stage1/stage1_sieve_gpu.c

NFS_GPU_OBJS = $(NFS_GPU_SRCS:.c=.no)

NFS_NOGPU_SRCS = \
	gnfs/poly/stage1/stage1_sieve_cpu.c

NFS_NOGPU_OBJS = $(NFS_NOGPU_SRCS:.c=.no)

ifeq ($(CUDA),1)
	NFS_HDR += $(NFS_GPU_HDR)
	NFS_SRCS += $(NFS_GPU_SRCS)
	NFS_OBJS += $(NFS_GPU_OBJS)
else
	NFS_HDR += $(NFS_NOGPU_HDR)
	NFS_SRCS += $(NFS_NOGPU_SRCS)
	NFS_OBJS += $(NFS_NOGPU_OBJS)
	GPU_OBJS =
endif

#---------------------------------- make targets -------------------------

help:
	@echo "to build:"
	@echo "make all"
	@echo "add 'WIN=1 if building on windows"
	@echo "add 'WIN64=1 if building on 64-bit windows"
	@echo "add 'ECM=1' if GMP-ECM is available (enables ECM)"
	@echo "add 'CUDA=1' for Nvidia graphics card support"
	@echo "add 'MPI=1' for parallel processing using MPI"
	@echo "add 'BOINC=1' to add BOINC wrapper"
	@echo "add 'NO_ZLIB=1' if you don't have zlib"

all: $(COMMON_OBJS) $(QS_OBJS) $(NFS_OBJS) $(GPU_OBJS)
	rm -f libmsieve.a
	ar r libmsieve.a $(COMMON_OBJS) $(QS_OBJS) $(NFS_OBJS)
	ranlib libmsieve.a
	$(CC) $(CFLAGS) demo.c -o msieve $(LDFLAGS) \
			libmsieve.a $(LIBS)

clean:
	cd cub && make clean WIN=$(WIN) WIN64=$(WIN64) && cd ..
	rm -f msieve msieve.exe libmsieve.a $(COMMON_OBJS) $(QS_OBJS) \
		$(NFS_OBJS) $(NFS_GPU_OBJS) $(NFS_NOGPU_OBJS) *.ptx

#----------------------------------------- build rules ----------------------

# common file build rules

%.o: %.c $(COMMON_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

# QS build rules

mpqs/sieve_core_generic_32k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -DBLOCK_KB=32 -DHAS_SSE2 \
		-DROUTINE_NAME=qs_core_sieve_generic_32k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_generic_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -DBLOCK_KB=64 -DHAS_SSE2 \
		-DROUTINE_NAME=qs_core_sieve_generic_64k \
		-c -o $@ mpqs/sieve_core.c

%.qo: %.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

# NFS build rules

%.no: %.c $(COMMON_HDR) $(NFS_HDR)
	$(CC) $(CFLAGS) -Ignfs -c -o $@ $<

# GPU build rules

stage1_core_sm20.ptx: $(NFS_GPU_HDR)
	$(NVCC) -arch sm_20 -ptx -o $@ $<

stage1_core_sm30.ptx: $(NFS_GPU_HDR)
	$(NVCC) -arch sm_30 -ptx -o $@ $<

stage1_core_sm35.ptx: $(NFS_GPU_HDR)
	$(NVCC) -arch sm_35 -ptx -o $@ $<

cub/built:
	cd cub && make WIN=$(WIN) WIN64=$(WIN64) sm=200,300,350 && cd ..
