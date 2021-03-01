CXX     = g++
#CXX		 = icpc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -Lhtslib -lz -lm -lbz2 -llzma -lpthread -fopenmp
INC = -Ihtslib
CXXFLAGS   = -g -Wall -O2 -fopenmp

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
INCLUDES+=	-I${VTUNE_HOME}/include
LIBS+=-L${VTUNE_HOME}/lib64 -littnotify
endif

all: debruijn.cpp common.cpp common.h
	$(CXX) debruijn.cpp common.cpp htslib/libhts.a $(CXXFLAGS) $(INC) $(LIBS) -o dbg

clean: 
	rm -f dbg
