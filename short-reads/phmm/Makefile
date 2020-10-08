CC      = gcc
CFLAGS= -g -O2 -std=c++11
GKLPATH= ../../tools/GKL
SHARED_LIBRARIES = -fopenmp -L$(GKLPATH)/build/native -lgkl_pairhmm_c -lstdc++

all: pairhmm

pairhmm: PairHMMUnitTest.o
	$(CC) -o $@ $(SHARED_LIBRARIES) $^ $(CFLAGS)

PairHMMUnitTest.o: PairHMMUnitTest.cpp PairHMMUnitTest.h
	$(CC) -c $(CFLAGS) -fopenmp $<

.PHONY: clean

clean:
	rm pairhmm *.o
