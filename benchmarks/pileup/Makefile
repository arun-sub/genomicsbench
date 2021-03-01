CC     = gcc
#CC		= icc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -fopenmp
CCFLAGS   = -g -Wall -pthread -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -O3 -fopenmp
INC = -I../../tools/htslib

ifneq ($(VTUNE_HOME),)
CCFLAGS += -DVTUNE_ANALYSIS=1
LIBS += -L${VTUNE_HOME}/lib64 -littnotity
INC += -I${VTUNE_HOME}/include
endif

all: medaka_common.c medaka_counts.c medaka_bamiter.c
	cd ../../tools/htslib && autoreconf -i && ./configure && $(MAKE)
	$(CC) $(CCFLAGS) medaka_common.c medaka_counts.c medaka_bamiter.c ../../tools/htslib/libhts.a $(INC) $(LIBS) -o pileup -std=c99 -O3

.PHONY: clean

clean: 
	cd ../../tools/htslib && $(MAKE) clean
	rm -f pileup

medaka_common.c: medaka_common.h
medaka_counts.c: medaka_counts.h medaka_common.h medaka_bamiter.h
medaka_bamiter.c: medaka_bamiter.h
