
all: inverse epsilon loop sum_numbers

inverse: inverse.f90
	gfortran -o $@ $^

epsilon: epsilon.f90
	gfortran -o $@ $^

loop: loop.f90
	gfortran -o $@ $^

sum_numbers: sum_numbers.f90
	gfortran -o $@ $^

clean:
	rm -f *~ inverse epsilon loop sum_numbers
