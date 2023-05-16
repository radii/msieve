
WIN = 0
WIN64 = 0

ifeq ($(WIN),1)
	CUDA_ROOT = $(shell echo $$CUDA_PATH)
	NVCC = "$(CUDA_ROOT)/bin/nvcc"
	EXT = dll
ifeq ($(WIN64),1)
	NVCCFLAGS += -m64
else
	NVCCFLAGS += -m32
endif
	NVCCFLAGS += -Xptxas -v -Xcudafe -\# -shared -Xptxas -abi=no
else
	NVCC = "$(shell which nvcc)"
	CUDA_ROOT = $(shell dirname $(NVCC))/../
	EXT = so
	NVCCFLAGS += -Xptxas -v -Xcudafe -\# -shared -Xptxas -abi=no \
			-Xcompiler -fPIC -Xcompiler -fvisibility=hidden
endif

COMMA = ,
ifdef sm
	SM_ARCH = $(subst $(COMMA),-,$(sm))
else 
    SM_ARCH = 200
endif

ifeq (520, $(findstring 520, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_52,code=\"sm_52,compute_52\" 
    SM_DEF 		+= -DSM520
endif
ifeq (370, $(findstring 370, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_37,code=\"sm_37,compute_37\" 
    SM_DEF 		+= -DSM370
endif
ifeq (350, $(findstring 350, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_35,code=\"sm_35,compute_35\" 
    SM_DEF 		+= -DSM350
endif
ifeq (300, $(findstring 300, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_30,code=\"sm_30,compute_30\"
    SM_DEF 		+= -DSM300
endif
ifeq (210, $(findstring 210, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_20,code=\"sm_21,compute_20\"
    SM_DEF 		+= -DSM210
endif
ifeq (200, $(findstring 200, $(SM_ARCH)))
    SM_TARGETS 	+= -gencode=arch=compute_20,code=\"sm_20,compute_20\"
    SM_DEF 		+= -DSM200
endif

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

CUB_DEPS = 	$(call rwildcard, cub,*.cuh)

INC = -I"$(CUDA_ROOT)/include" -I.

DEPS = ./Makefile \
	sort_engine.cu \
	sort_engine.h \
	$(CUB_DEPS)

LIBNAME = sort_engine

all: $(LIBNAME).$(EXT)
	touch built

clean :
	rm -f  *.dll *.so *.lib *.exp built

$(LIBNAME).$(EXT) : $(DEPS)
	$(NVCC) $(SM_TARGETS) $(SM_DEF) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

