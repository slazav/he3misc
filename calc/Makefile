#FFLAGS= -Werror -Wconversion\ 
#  -Wintrinsic-shadow -Wline-truncation\ 
#  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\ 
#  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\ 
#  -Wno-unused-parameter -Wno-align-commons 

# gfortran parameters
FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter

TARGETS=he3 he3dat
all: $(TARGETS)

FC=gfortran

clean:
	rm -f $(TARGETS) *.o *.a

he3: he3.f ../lib/libhe3.a
vft: vft.f


CFLAGS=-g -I ../lib
FFLAGS=-g -I ../lib

he3dat: he3dat.c he3dat.h
	$(CC) he3dat.c ../lib/libhe3.a -lgfortran -lc -lm -o he3dat
#he3_funcs.o: he3_funcs.f
