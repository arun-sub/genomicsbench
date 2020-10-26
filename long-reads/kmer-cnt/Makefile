#CXX     = g++
CXX		 = icpc
VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -lz -lm -lbz2 -llzma -L${VTUNE_HOME}/lib64 -littnotify -fopenmp
CXXFLAGS   = -O3 -fopenmp
INC = -Ilibcuckoo -I${VTUNE_HOME}/include


all: sequence_container.cpp sequence.cpp vertex_index.cpp kmer_cnt.cpp
	$(CXX) $(CXXFLAGS) sequence_container.cpp sequence.cpp vertex_index.cpp kmer_cnt.cpp $(INC) $(LIBS) -o kmer-cnt

.PHONY: clean

clean: 
	rm -f kmer_cnt

sequence_container.cpp: sequence_container.h sequence.h logger.h
sequence.cpp: sequence.h
vertex_index.cpp: vertex_index.h parallel.h memory_info.h logger.h config.h sequence_container.h kmer.h
kmer-cnt.cpp: sequence_container.h sequence.h vertex_index.h memory_info.h logger.h utils.h config.h parallel.h
