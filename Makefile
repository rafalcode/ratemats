# Programs for subsitution rate matrices
#
CC=gcc
CFLAGS=-g -Wall# -pg # note the gprof option
EXES=ratematn
SPECLIBS=-lcairo
GSL_LIBS=-lgsl -lgslcblas -lm

# in a previous incarnation, this was the dric7.c program
# the n at the end implies numbers. It's to check validity
# so this will generate no graphics
ratematn: ratematn.c
	@${CC} ${CFLAGS} -o $@ $^ -DUNPREDRA -lm
ratematn_d: ratematn.c
	@${CC} ${CFLAGS} -DDBG -o $@ $^ -DUNPREDRA -lm

.PHONY: clean

clean:
	rm -f ${EXES}
