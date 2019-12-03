SHELL=/bin/sh

all: do-hello do-exp do-visible

clean:
	rm -f *.s *.pp.c *.o hello exp *.mod exp-?

do-hello: hello.pp.c hello.s hello.o hello-s.o hello

hello.pp.c: hello.c
	gcc -E -o $@ $<

hello.s: hello.c
	gcc -S $<

hello.o: hello.c
	gcc -c $<

hello-s.o: hello.s
	as -o $@ $<

hello: hello.o
	gcc -o $@ $<


do-exp: exp.pp.c exp.s exp.o exp-s.o exp-g.o exp exp-m

exp.pp.c: exp.c
	gcc -E -o $@ $<

exp.s: exp.c
	gcc -S $<

exp.o: exp.c
	gcc -c $<

exp-g.o: exp.c
	gcc -g -o $@ -c $<

exp-s.o: exp.s
	as -o $@ $<

exp: exp.o
	-gcc -o $@ $<

exp-m: exp.o
	-gcc -o $@ $< -lm


do-visible: visibility-c.o visibility-cpp.o visibility-f90.o

visibility-c.o: visibility.c
	gcc -c $< -o $@

visibility-cpp.o: visibility.cpp
	g++ -c $< -o $@

visibility-f90.o: visibility.f90
	gfortran -c $< -o $@


