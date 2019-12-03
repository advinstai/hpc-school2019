FC       = gfortran 
#FC       = ifort
#FC        = pgf77

CFLAGS	= $(ABI) $(ARCH)

TARGETS	= prog1 prog2a prog2b prog2c prog2d prog3

default: all

prog1: prog1.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O3 prog1.f

prog2a: prog2.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O3 prog2.f

prog2b: prog2.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O0 prog2.f

prog2c: prog2.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O3 -ffast-math prog2.f

prog2d: prog2.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O3 -ffast-math -funsafe-math-optimizations prog2.f

prog3: prog3.f
	$(FC) -o $@ $(ABI) $(ARCH) \
	      -O3 prog3.f
all: $(TARGETS)

clean:
	/bin/rm -f $(TARGETS)
	/bin/rm -f *.s *.w2f.f
