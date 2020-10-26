#CXX=g++
CXX=icpc
CXXFLAGS=-O3 -fopenmp
VTUNE_HOME= /opt/intel/vtune_profiler
INCLUDES=-I../../tools/spoa/include/ -I${VTUNE_HOME}/include
LIBS=-L../../tools/spoa/build/lib/ -lspoa -L${VTUNE_HOME}/lib64 -littnotify -fopenmp

all: msa_spoa_omp.cpp
	$(CXX) $(CXXFLAGS) msa_spoa_omp.cpp ${INCLUDES} ${LIBS} -o msa_spoa_omp

.PHONY: clean

clean:
	rm -f msa_spoa_omp
