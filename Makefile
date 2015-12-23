CC=gcc
EXEEXT=
#EXEEXT=.exe     #for windows
CFLAGS=-Wall
gmmtest$(EXEEXT): gmm.o gmmtest.o
	$(CC) -o $@ $^ -lm

gmm.o: gmm.c gmm.h
gmmtest.o: gmmtest.c gmm.h

