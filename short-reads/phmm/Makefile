#CC      = gcc
CC		= icpc
CFLAGS= -g -O2 -std=c++11
VTUNE_HOME= /opt/intel/vtune_profiler
INC		= -I${VTUNE_HOME}/include
GKLPATH= ../../tools/GKL
SHARED_LIBRARIES = -fopenmp -L$(GKLPATH)/build/native -lgkl_pairhmm_c -L${VTUNE_HOME}/lib64 -littnotify -lstdc++

all: pairhmm

pairhmm: PairHMMUnitTest.o
	$(CC) -o $@ $^ $(CFLAGS) ${INC} $(SHARED_LIBRARIES)

PairHMMUnitTest.o: PairHMMUnitTest.cpp PairHMMUnitTest.h
	$(CC) -c $(CFLAGS) ${INC} -fopenmp $<

.PHONY: clean

clean:
	rm -f pairhmm *.o
