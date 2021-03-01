CXX=g++
#CXX=icpc
CXXFLAGS=-O3 -fopenmp
#VTUNE_HOME= /opt/intel/vtune_profiler
INCLUDES=-I../../tools/spoa/include
LIBS=-L../../tools/spoa/build/lib/ -lspoa -fopenmp

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
