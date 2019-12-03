CC       = gcc

CFLAGS	= $(ABI) $(ARCH)

TARGETS	= prog1 prog2 prog3 prog4 prog5
SRC = paranoia.c args.h

default: all

prog1:  $(SRC)
	$(CC) -o $@ $(ABI) $(ARCH) \
	      -O0 -mfpmath=387 paranoia.c -lm

prog2:  $(SRC)
	$(CC) -o $@ $(ABI) $(ARCH) \
	      -O0 -mfpmath=sse paranoia.c -lm

prog3:  $(SRC)
	$(CC) -o $@ $(ABI) $(ARCH) \
	      -O0 -mpc80 paranoia.c -lm

prog4:  $(SRC)
	$(CC) -o $@ $(ABI) $(ARCH) \
	      -O2 -mpc64 -ansi -pedantic paranoia.c -lm

prog5:  $(SRC)
	$(CC) -o $@ $(ABI) $(ARCH) \
	      -O3 -ffast-math -funsafe-math-optimizations paranoia.c -lm

all: $(TARGETS)

clean:
	/bin/rm -f $(TARGETS)
	/bin/rm -f *.s *.w2f.f
