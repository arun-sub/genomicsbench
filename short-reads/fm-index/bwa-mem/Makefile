BWA_PATH=../../../tools/bwa

EXE=fmi_smem_reseed_test
CC=gcc
CFLAGS=-g -Wall -Wextra -O2 -fopenmp
DFLAGS=
INCLUDES=-I${BWA_PATH}
LIBS=-L${BWA_PATH} -fopenmp -lm -lpthread -lz -lbwa


.PHONY:all clean depend
.SUFFIXES:.c .o


.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

fmi_smem_reseed_test: ${BWA_PATH}/libbwa.a fmi_smem_reseed_test.o
		$(CC) -o $@ $^ $(LIBS)

clean:
		rm -fr *.o $(EXE)

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

fmi_smem_reseed_test.o: ${BWA_PATH}/bntseq.h ${BWA_PATH}/bwt.h ${BWA_PATH}/utils.h ${BWA_PATH}/bwa.h
fmi_smem_reseed_test.o: ${BWA_PATH}/kseq.h 
