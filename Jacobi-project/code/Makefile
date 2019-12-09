# inputs for the executables

dim = 10
it  = 100
r = 3
c = 3

CC=cc
CFLAGS=-O3

SRCS=$(wildcard *.c)
EXE=$(SRCS:.c=.x)

all: $(EXE)

run: clean all
	./$(EXE) $(dim) $(it) $(r) $(c)

%.x: %.c
	$(CC) $< $(CFLAGS) -o $@

clean:
	@rm -f *~ $(EXE) solution.dat 

plot:
	@gnuplot -p plot.plt

.PHONY: clean plot all

