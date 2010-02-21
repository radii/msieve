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

# gcc with basic optimization (-march flag could
# get overridden by architecture-specific builds)
CC = gcc -D_FILE_OFFSET_BITS=64
WARN_FLAGS = -Wall -W
OPT_FLAGS = -O3 -fomit-frame-pointer -march=athlon-xp -DNDEBUG
#OPT_FLAGS = -O3 -fomit-frame-pointer -march=k8 -DNDEBUG

CFLAGS = $(OPT_FLAGS) $(MACHINE_FLAGS) $(WARN_FLAGS) \
		-I. -Iinclude -Ignfs -Ignfs/poly -Ignfs/poly/stage1

# tweak the compile flags

ifeq ($(ECM),1)
	CFLAGS += -DHAVE_GMP_ECM
	LIBS += -lecm
endif
ifeq ($(WIN32_3GB),1)
	LDFLAGS += -Wl,--large-address-aware
endif
ifeq ($(CUDA),1)
	# these environment variables are set in windows
	# but not in linux; attempt to sample them anyway
	CUDA_INC_DIR = $(shell echo $$CUDA_INC_PATH)
	CUDA_LIB_DIR = $(shell echo $$CUDA_LIB_PATH)
	CFLAGS += -I"$(CUDA_INC_DIR)" -DHAVE_CUDA

	# Also, the CUDA driver library has a different name in linux
	LIBS += "$(CUDA_LIB_DIR)/cuda.lib"
	# LIBS += -lcuda
endif

# Note to MinGW users: the library does not use pthreads calls in
# win32 or win64, so it's safe to pull libpthread into the link line.
# Of course this does mean you have to install the minGW pthreads bundle...

LIBS += -lgmp -lm -lpthread

#---------------------------------- Generic file lists -------------------

COMMON_HDR = \
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
	include/util.h

COMMON_SRCS = \
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
	common/mp.c \
	common/polyroot.c \
	common/prime_delta.c \
	common/prime_sieve.c \
	common/savefile.c \
	common/strtoll.c \
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
	mpqs/sqrt.qo

QS_CORE_OBJS = \
	mpqs/sieve_core_generic_32k.qo \
	mpqs/sieve_core_generic_64k.qo

QS_CORE_OBJS_X86 = \
	mpqs/sieve_core_p2_64k.qo \
	mpqs/sieve_core_p3_64k.qo \
	mpqs/sieve_core_p4_64k.qo \
	mpqs/sieve_core_pm_32k.qo \
	mpqs/sieve_core_core_32k.qo \
	mpqs/sieve_core_k7_64k.qo \
	mpqs/sieve_core_k7xp_64k.qo \
	mpqs/sieve_core_k8_64k.qo

QS_CORE_OBJS_X86_64 = \
	mpqs/sieve_core_p4_64_64k.qo \
	mpqs/sieve_core_core_64_32k.qo \
	mpqs/sieve_core_k8_64_64k.qo

#---------------------------------- GPU file lists -------------------------

GPU_OBJS = \
	stage1_core_deg46_48.ptx \
	stage1_core_deg46_64.ptx \
	stage1_core_deg5_128.ptx \
	stage1_core_deg5_48.ptx \
	stage1_core_deg5_64.ptx \
	stage1_core_deg5_72.ptx \
	stage1_core_deg5_96.ptx \
	stage1_core_deg6_128.ptx \
	stage1_core_deg6_72.ptx \
	stage1_core_deg6_96.ptx

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
	gnfs/poly/stage1/stage1_core_gpu/cuda_intrinsics.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg46_64.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg5_128.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg5_64.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg5_96.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg6_96.h \
	gnfs/poly/stage1/stage1_core_gpu/stage1_core_deg6_128.h

NFS_NOGPU_HDR = \
	gnfs/poly/stage1/stage1_core/cpu_intrinsics.h

NFS_SRCS = \
	gnfs/poly/poly.c \
	gnfs/poly/poly_skew.c \
	gnfs/poly/polyutil.c \
	gnfs/poly/root_score.c \
	gnfs/poly/size_score.c \
	gnfs/poly/stage1/stage1.c \
	gnfs/poly/stage1/stage1_roots.c \
	gnfs/poly/stage1/stage1_sieve.c \
	gnfs/poly/stage2/optimize.c \
	gnfs/poly/stage2/root_sieve.c \
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
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg46_64.c \
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg5_128.c \
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg5_64.c \
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg5_96.c \
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg6_128.c \
	gnfs/poly/stage1/stage1_core_gpu/stage1_sieve_deg6_96.c

NFS_GPU_OBJS = $(NFS_GPU_SRCS:.c=.no)

NFS_NOGPU_SRCS = \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg46_64.c \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg5_128.c \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg5_64.c \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg5_96.c \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg6_128.c \
	gnfs/poly/stage1/stage1_core/stage1_sieve_deg6_96.c

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

all:
	@echo "pick a target:"
	@echo "x86       32-bit Intel/AMD systems (required if gcc used)"
	@echo "x86_64    64-bit Intel/AMD systems (required if gcc used)"
	@echo "generic   portable code"
	@echo "add 'ECM=1' if GMP-ECM is available (enables ECM)"
	@echo "add 'CUDA=1' for Nvidia graphics card support"

x86: $(COMMON_OBJS) $(QS_OBJS) $(QS_CORE_OBJS) \
		$(QS_CORE_OBJS_X86) $(NFS_OBJS) $(GPU_OBJS)
	rm -f libmsieve.a
	ar r libmsieve.a $(COMMON_OBJS) $(QS_OBJS) \
			$(QS_CORE_OBJS) $(QS_CORE_OBJS_X86) \
			$(NFS_OBJS)
	ranlib libmsieve.a
	$(CC) $(CFLAGS) demo.c -o msieve $(LDFLAGS) \
			libmsieve.a $(LIBS)

x86_64: $(COMMON_OBJS) $(QS_OBJS) $(QS_CORE_OBJS) \
		$(QS_CORE_OBJS_X86_64) $(NFS_OBJS) $(GPU_OBJS)
	rm -f libmsieve.a
	ar r libmsieve.a $(COMMON_OBJS) $(QS_OBJS) \
			$(QS_CORE_OBJS) $(QS_CORE_OBJS_X86_64) \
			$(NFS_OBJS)
	ranlib libmsieve.a
	$(CC) $(CFLAGS) demo.c -o msieve $(LDFLAGS) \
			libmsieve.a $(LIBS)

generic: $(COMMON_OBJS) $(QS_OBJS) $(QS_CORE_OBJS) $(NFS_OBJS) $(GPU_OBJS)
	rm -f libmsieve.a
	ar r libmsieve.a $(COMMON_OBJS) $(QS_OBJS) \
			$(QS_CORE_OBJS) $(NFS_OBJS)
	ranlib libmsieve.a
	$(CC) $(CFLAGS) demo.c -o msieve $(LDFLAGS) \
			libmsieve.a $(LIBS)

clean:
	rm -f msieve msieve.exe libmsieve.a $(COMMON_OBJS) 	\
		$(QS_OBJS) $(QS_CORE_OBJS) $(QS_CORE_OBJS_X86) \
		$(QS_CORE_OBJS_X86_64) $(NFS_OBJS) \
		$(NFS_GPU_OBJS) $(NFS_NOGPU_OBJS) *.ptx

#----------------------------------------- build rules ----------------------

# common file build rules

%.o: %.c $(COMMON_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

# QS build rules

mpqs/sieve_core_generic_32k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -DBLOCK_KB=32 -DCPU_GENERIC \
		-DROUTINE_NAME=qs_core_sieve_generic_32k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_generic_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -DBLOCK_KB=64 -DCPU_GENERIC \
		-DROUTINE_NAME=qs_core_sieve_generic_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_p2_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=pentium2 -DBLOCK_KB=64 -DCPU_PENTIUM2 \
		-DROUTINE_NAME=qs_core_sieve_p2_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_p3_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=pentium3 -DBLOCK_KB=64 -DCPU_PENTIUM3 \
		-DROUTINE_NAME=qs_core_sieve_p3_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_p4_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=pentium4 -DBLOCK_KB=64 -DCPU_PENTIUM4 \
		-DROUTINE_NAME=qs_core_sieve_p4_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_pm_32k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=pentium-m -DBLOCK_KB=32 -DCPU_PENTIUM_M \
		-DROUTINE_NAME=qs_core_sieve_pm_32k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_core_32k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=prescott -DBLOCK_KB=32 -DCPU_CORE \
		-DROUTINE_NAME=qs_core_sieve_core_32k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_k7_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=athlon -DBLOCK_KB=64 -DCPU_ATHLON \
		-DROUTINE_NAME=qs_core_sieve_k7_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_k7xp_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=athlon-xp -DBLOCK_KB=64 -DCPU_ATHLON_XP \
		-DROUTINE_NAME=qs_core_sieve_k7xp_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_k8_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=k8 -DBLOCK_KB=64 -DCPU_OPTERON \
		-DROUTINE_NAME=qs_core_sieve_k8_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_p4_64_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=nocona -DBLOCK_KB=64 -DCPU_PENTIUM4 \
		-DROUTINE_NAME=qs_core_sieve_p4_64k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_core_64_32k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=nocona -DBLOCK_KB=32 -DCPU_CORE \
		-DROUTINE_NAME=qs_core_sieve_core_32k \
		-c -o $@ mpqs/sieve_core.c

mpqs/sieve_core_k8_64_64k.qo: mpqs/sieve_core.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -march=k8 -DBLOCK_KB=64 -DCPU_OPTERON \
		-DROUTINE_NAME=qs_core_sieve_k8_64k \
		-c -o $@ mpqs/sieve_core.c

%.qo: %.c $(COMMON_HDR) $(QS_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

# NFS build rules

%.no: %.c $(COMMON_HDR) $(NFS_HDR)
	$(CC) $(CFLAGS) -Ignfs -c -o $@ $<

# GPU build rules

%.ptx: gnfs/poly/stage1/stage1_core_gpu/%.cu $(NFS_GPU_HDR)
	nvcc -ptx -o $@ $<
