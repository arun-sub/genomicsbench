CC      = gcc
#CC		= icpc
CFLAGS= -g -O2 -std=c++11
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
