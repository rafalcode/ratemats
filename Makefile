# Programs for subsitution rate matrices
# Mne: RATEMAT, an "n" after means numbers only.
#
CC=gcc
CFLAGS=-g -Wall# -pg # note the gprof option
EXES=ratematn ratematn_d bratematn bratematn_d discrate discrate_d waitimes
SPECLIBS=-lcairo
GSL_LIBS=-lgsl -lgslcblas -lm

# in a previous incarnation, this was the dric7.c program
# the n at the end implies numbers. It's to check validity
# so this will generate no graphics
ratematn: ratematn.c
	@${CC} ${CFLAGS} -o $@ $^ -lm
ratematn_d: ratematn.c
	@${CC} ${CFLAGS} -DDBG -o $@ $^ -lm

# binary case, only two symbols. CONTINUOUS case
	@${CC} ${CFLAGS} -o $@ $^ -lm
bratematn_d: bratematn.c
	@${CC} ${CFLAGS} -DDBG -o $@ $^ -lm

# binary, and, discrete rate case, only two symbols ... and, DISCRETE VERSION
bdiscrate: bdiscrate.c
	@${CC} ${CFLAGS} -o $@ $^ -lm
bdiscrate_d: bdiscrate.c
	@${CC} ${CFLAGS} -DDBG -o $@ $^ -lm

# proto version
# actually for comparing different "wrong" versions of the exponential waiting time
ratematproto: ratematproto.c
	@${CC} ${CFLAGS} -o $@ $^ -lm
ratematproto_d: ratematproto.c
	@${CC} ${CFLAGS} -DDBG -o $@ $^ -lm

waitimes: waitimes.c
	@${CC} ${CFLAGS} -o $@ $^ -lm

.PHONY: clean

clean:
	rm -f ${EXES}
