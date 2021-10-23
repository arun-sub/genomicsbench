##/*************************************************************************************
##                           The MIT License
##
##   Copyright (C) 2020 Intel Labs.
##
##   Permission is hereby granted, free of charge, to any person obtaining a copy
##   of this software and associated documentation files (the "Software"), to deal
##   in the Software without restriction, including without limitation the rights
##   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##   copies of the Software, and to permit persons to whom the Software is
##   furnished to do so, subject to the following conditions:
##
##   The above copyright notice and this permission notice shall be included in all
##   copies or substantial portions of the Software.
##   
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>
##		Saurabh Kalikar <saurabh.kalikar@intel.com>
##*****************************************************************************************/

BENCH_SMEM=		smem
#CXX=		icpc
ifeq ($(CXX), icpc)
    CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
endif

ARCH_FLAGS=	-msse4.1
#ARCH_FLAGS=	-mavx512bw
MEM_FLAGS=	-DSAIS=1
CPPFLAGS=	-DENABLE_PREFETCH $(MEM_FLAGS) -DKSW=1
INCLUDES=   -Iext -Iext/safestringlib/include -I.
LIBS=		-fopenmp -lm -lz -L. -ltal -Lext/safestringlib/ -lsafestring
OBJS=		ext/utils.o \
			ext/kstring.o ext/bntseq.o \
			FMI_search.o ext/bseq.o
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a
TAL_LIB=    libtal.a

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

#CXXFLAGS=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2
CXXFLAGS=	-g -O3 -fopenmp $(ARCH_FLAGS)

#.PHONY:all clean depend
.PHONY:all clean depend
#.PHONY:all clean
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(BENCH_SMEM)

$(BENCH_SMEM):$(TAL_LIB) $(SAFE_STR_LIB) bench-smem.o
	$(CXX) $(CXXFLAGS) bench-smem.o $(LIBS) -o $@

$(TAL_LIB):$(OBJS)
		ar rcs $(TAL_LIB) $(OBJS)

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && make clean && make CC=$(CC) directories libsafestring.a

clean:
	rm -fr *.o ext/*.o $(BENCH_SMEM) $(TAL_LIB)
	cd ext/safestringlib && make CC=$(CC) clean

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)

# DO NOT DELETE

FMI_search.o: FMI_search.h ext/bntseq.h
FMI_search.o: ext/utils.h  ext/sais.h
ext/bntseq.o: ext/bntseq.h ext/utils.h ext/kseq.h
#ext/bwa.o: ext/bntseq.h ext/bwa.h ext/utils.h
#ext/bwa.o: ext/kstring.h ext/kseq.h
ext/kstring.o: ext/kstring.h
ext/utils.o: ext/utils.h ext/kseq.h
ext/bseq.o: ext/bseq.h ext/utils.h ext/kseq.h
benchmarks/bench-smem.o: FMI_search.h ext/bntseq.h ext/utils.h
benchmarks/bench-smem.o:  ext/sais.h
