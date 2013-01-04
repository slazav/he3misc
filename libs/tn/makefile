all: tn.a

GOBJS =	blas.o tn.o
FC = gfortran

.SUFFIXES: .f .o

.f.o:
	gfortran -c $<

tn.a: $(GOBJS)
	ar rs $@ $+

clean:
	rm -f -- *.o *.a
