CC     = gcc
LIBS = -Lhtslib -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto
CCFLAGS   = -g -Wall -pthread -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -O3
INC = -Ihtslib

all: medaka_common.c medaka_counts.c medaka_bamiter.c
	$(CC) $(CCFLAGS) medaka_common.c medaka_counts.c medaka_bamiter.c htslib/libhts.a $(INC) $(LIBS) -o pileup -std=c99 -O3

.PHONY: clean

clean: 
	rm -f pileup

medaka_common.c: medaka_common.h
medaka_counts.c: medaka_counts.h medaka_common.h medaka_bamiter.h
medaka_bamiter.c: medaka_bamiter.h
