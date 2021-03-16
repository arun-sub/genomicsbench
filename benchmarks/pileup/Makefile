CC     = gcc
#CC		= icc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ldl
CFLAGS   = -g -Wall -pthread -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -O3 -fopenmp
INC = -I../../tools/htslib

ifneq ($(VTUNE_HOME),)
CFLAGS += -DVTUNE_ANALYSIS=1
LIBS += -L${VTUNE_HOME}/lib64 -littnotify
INC += -I${VTUNE_HOME}/include
endif

all: medaka_common.c medaka_counts.c medaka_bamiter.c
	$(CC) medaka_common.c medaka_counts.c medaka_bamiter.c ../../tools/htslib/libhts.a $(CFLAGS) $(INC) $(LIBS) -o pileup

.PHONY: clean

clean: 
	#cd ../../tools/htslib && $(MAKE) clean
	rm -f pileup

medaka_common.c: medaka_common.h
medaka_counts.c: medaka_counts.h medaka_common.h medaka_bamiter.h
medaka_bamiter.c: medaka_bamiter.h
