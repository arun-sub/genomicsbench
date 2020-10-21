CXX=g++
CXXFLAGS=-O3 -fopenmp
INCLUDES=../../tools/spoa/include/
LIBS=../../tools/spoa/build/lib/ -lspoa -fopenmp

all: msa_spoa_omp.cpp
	$(CXX) $(CXXFLAGS) msa_spoa_omp.cpp -I${INCLUDES} -L${LIBS} -o msa_spoa_omp

.PHONY: clean

clean:
	rm -f msa_spoa_omp
