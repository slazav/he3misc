# Makefile

# SHELL=/bin/sh

# Fortran compiler, compiler flags, linker flags
# object files and name of executable

FC=gfortran
OBJS= mainlib.o general.o free.o texture.o modu.o gap.o profiles.o nmr.o
OBJ2=mainprog.o
EXE=texture
LIBR=texture.so
FFLAGS=-O
FFLAGX=-fPIC

USEBTN=0

ifeq ($(USEBTN),1)
  TNLIB=../libs/711/btn.a
else
  TNLIB=../libs/tn/tn.a
endif


all:	$(EXE) $(LIBR)

$(EXE):		$(OBJ2) $(LIBR) $(TNLIB)
		$(FC) -g -o $(EXE) $+ $(LDFLAGS) -Wl,-rpath=$(PWD)

mainprog.o:	mainprog.f90
		$(FC) -g -c mainprog.f90 $(FFLAGS)

$(LIBR):	$(OBJS)
		$(FC) -g -shared -o $(LIBR) $(OBJS) $(LDFLAGS) $(FFLAGS)


mainlib.o:	mainlib.f90 texture.o general.o modu.o free.o gap.o profiles.o nmr.o
		$(FC) -g -cpp -DUSEBTN=$(USEBTN) -c mainlib.f90 $(FFLAGS) $(FFLAGX)


free.o:		free.f90 texture.o general.o modu.o 
		$(FC) -g -c free.f90 $(FFLAGS) $(FFLAGX)


nmr.mod nmr.o:	nmr.f90 general.mod
		$(FC) -g -c nmr.f90 $(FFLAGS) $(FFLAGX)


modu.mod modu.o: modu.f90 general.o
		$(FC) -g -c modu.f90 $(FFLAGS) $(FFLAGX)

profiles.o:	profiles.f90 general.o modu.o free.o
		$(FC) -g -c profiles.f90 $(FFLAGS) $(FFLAGX)


general.mod general.o:  general.f90
		$(FC) -g -c general.f90 $(FFLAGS) $(FFLAGX)

texture.o:	texture.f90 gap.o general.o
		$(FC) -g -c texture.f90 $(FFLAGS) $(FFLAGX)

gap.o:		gap.f90 general.o
		$(FC) -g -c gap.f90 $(FFLAGS) $(FFLAGX)

blas.o:		blas.f 
		$(FC) -g -c -FI -zero -fno-automatic -std=legacy blas.f $(FFLAGS) $(FFLAGX)

tn.o:		tn.f
		$(FC) -g -c -FI -zero -fno-automatic -std=legacy tn.f $(FFLAGS) $(FFLAGX)

clean:
	rm -f *.o *.mod *.so texture alpha.dat beta.dat spec.dat