CXX=g++
#CXX=icpc


ifeq ($(arch),sse41)
	ARCH_FLAGS=-msse4.1
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else	
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else	
		ARCH_FLAGS=-mavx512bw
	endif
else ifeq ($(arch),native)
		ARCH_FLAGS=-march=native
else ifneq ($(arch),)
	## To provide a different architecture flag like -march=core-avx2.
		ARCH_FLAGS=$(arch)
endif

CXXFLAGS=-O3 -fopenmp $(ARCH_FLAGS)
#VTUNE_HOME= /opt/intel/vtune_profiler
INCLUDES=-I../../tools/spoa/include
LIBS=-L../../tools/spoa/build/lib/ -lspoa -ldl -fopenmp

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
INCLUDES+=	-I${VTUNE_HOME}/include
LIBS+=-L${VTUNE_HOME}/lib64 -littnotify
endif

all: msa_spoa_omp.cpp
	$(CXX) $(CXXFLAGS) msa_spoa_omp.cpp ${INCLUDES} ${LIBS} -o poa

.PHONY: clean

clean:
	rm -f poa
