CC=gcc
CFLAGS=-g -Wall -Wextra -Wc++-compat -O2
CPPFLAGS=

INCLUDES=-I.

OBJS=ksw.o

PROG=ksw-test
LIBS=-lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(OBJS)
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

ksw.o: ksw.h kseq.h
