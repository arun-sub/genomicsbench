CXX     = g++
#CXX		 = icpc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -L../../tools/htslib -lcurl -lcrypto -lz -lm -lbz2 -llzma -lpthread -ldl -fopenmp
INC = -I../../tools/htslib
CXXFLAGS   = -g -Wall -O2 -fopenmp

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
INC+=	-I${VTUNE_HOME}/include
LIBS+=-L${VTUNE_HOME}/lib64 -littnotify
endif

all: debruijn.cpp common.cpp common.h
	$(CXX) debruijn.cpp common.cpp ../../tools/htslib/libhts.a $(CXXFLAGS) $(INC) $(LIBS) -o dbg

clean: 
	rm -f dbg
