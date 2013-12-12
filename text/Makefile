# Makefile

# SHELL=/bin/sh

# Fortran compiler, compiler flags, linker flags
# object files and name of executable

# gfortran parameters
FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O

FC=gfortran

LDFLAGS=-L../lib/external

OBJS= mainlib.o general.o free.o modu.o profiles.o nmr.o
OBJ2=mainprog.o
EXE=texture
LIBR=texture.so
FFLAGX=$(FFLAGS) -fPIC

USEBTN=0

ifeq ($(USEBTN),0)
  TNLIB=-ltn
else
  TNLIB=-lbtn
endif


all:	$(EXE) $(LIBR)

$(EXE):		$(OBJ2) $(LIBR) ../lib/he3_const.o ../lib/libhe3.a
		$(FC) -g -o $(EXE) $+ $(TNLIB) $(LDFLAGS) -Wl,-rpath=$(PWD) $(FFLAGS)

mainprog.o:	mainprog.c
		$(CC) -g -c mainprog.c

$(LIBR):	$(OBJS)
		$(FC) -g -shared -o $(LIBR) $(OBJS) $(LDFLAGS) $(FFLAGS)


mainlib.o:	mainlib.f90 general.o modu.o free.o profiles.o nmr.o
		$(FC) -g -cpp -DUSEBTN=$(USEBTN) -c mainlib.f90 $(FFLAGS)


free.o:		free.f90 general.o modu.o 
		$(FC) -g -c free.f90 $(FFLAGS)


nmr.mod nmr.o:	nmr.f90 general.mod
		$(FC) -g -c nmr.f90 $(FFLAGS)


modu.mod modu.o: modu.f90 general.o
		$(FC) -g -c modu.f90 $(FFLAGS)

profiles.o:	profiles.f90 general.o modu.o free.o
		$(FC) -g -c profiles.f90 $(FFLAGS)


general.mod general.o:  general.f90
		$(FC) -g -c general.f90 $(FFLAGS)

blas.o:		blas.f 
		$(FC) -g -c -FI -zero -fno-automatic -std=legacy blas.f $(FFLAGS)

tn.o:		tn.f
		$(FC) -g -c -FI -zero -fno-automatic -std=legacy tn.f $(FFLAGS)
	)

clean:
	rm -f *.o *.mod *.so texture alpha.dat beta.dat spec.dat