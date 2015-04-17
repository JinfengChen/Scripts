.SUFFIXES: .c .cpp

CC=g++
CFLAGS = -g 

%.o:	%.cpp
	$(CC) -c -o $@ $(CFLAGS) $<

all: vit trainhmm genseq

vit:	hmm.o tables.o logprobs.o testviterbi.o
	$(CC) -o $@ $(CFLAGS) hmm.o tables.o logprobs.o testviterbi.o -lm

trainhmm:	hmm.o tables.o logprobs.o trainhmm.o
	$(CC) -o $@ $(CFLAGS) hmm.o tables.o logprobs.o trainhmm.o -lm

genseq:	hmm.o tables.o logprobs.o genseq.o
	$(CC) -o $@ $(CFLAGS) hmm.o tables.o logprobs.o genseq.o -lm

clean:
	rm -f *.o

depend: 
	g++ -M *.cpp >.dep

include .dep
