CC      = gcc
#CC		= icpc
ARCH_FLAGS=	-msse4.1

ifeq ($(arch),sse)
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
CFLAGS= -g -O2 -std=c++11 $(ARCH_FLAGS)

#VTUNE_HOME= /opt/intel/vtune_profiler
GKLPATH= ../../tools/GKL
INC=
SHARED_LIBRARIES = -fopenmp -L$(GKLPATH)/build/native -lgkl_pairhmm_c -lstdc++

ifneq ($(VTUNE_HOME),)
CFLAGS+= -DVTUNE_ANALYSIS=1
INC+= -I${VTUNE_HOME}/include
SHARED_LIBRARIES+= -L${VTUNE_HOME}/lib64 -littnotify
endif

all: phmm

phmm: PairHMMUnitTest.o
	$(CC) -o $@ $^ $(CFLAGS) ${INC} $(SHARED_LIBRARIES)

PairHMMUnitTest.o: PairHMMUnitTest.cpp PairHMMUnitTest.h
	$(CC) -c $(CFLAGS) ${INC} -fopenmp $<

.PHONY: clean

clean:
	rm -f phmm *.o
