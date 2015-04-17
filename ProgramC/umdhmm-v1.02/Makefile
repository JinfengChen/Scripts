#
#
# Make file for compiling HMM code in this directory.
# Author: Tapas Kanungo
# Date: 23 February 1998
# $Id: Makefile,v 1.3 1998/02/23 08:12:35 kanungo Exp kanungo $
# 
#
CFLAGS= -g
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
SRCS=baum.c viterbi.c forward.c backward.c hmmutils.c sequence.c \
	genseq.c nrutil.c testvit.c esthmm.c hmmrand.c testfor.c 

all :	genseq testvit testfor esthmm
	
genseq: genseq.o sequence.o nrutil.o hmmutils.o  hmmrand.o
	 $(CC) -o genseq genseq.o sequence.o nrutil.o \
	hmmrand.o hmmutils.o  -lm
testvit: testvit.o viterbi.o nrutil.o hmmutils.o sequence.o 
	 $(CC) -o testvit testvit.o viterbi.o nrutil.o sequence.o \
		hmmutils.o  hmmrand.o -lm
testfor: testfor.o forward.o nrutil.o hmmutils.o sequence.o hmmrand.o
	 $(CC) -o testfor testfor.o forward.o nrutil.o sequence.o \
		hmmutils.o  hmmrand.o -lm
esthmm: esthmm.o baum.o nrutil.o hmmutils.o sequence.o \
		forward.o backward.o hmmrand.o
	 $(CC) -o esthmm esthmm.o baum.o nrutil.o sequence.o hmmutils.o \
		forward.o backward.o hmmrand.o -lm
clean:
	rm *.o a.out 
# DO NOT DELETE THIS LINE -- make depend depends on it.

