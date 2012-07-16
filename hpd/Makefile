FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

TARGETS=vmcw vft he3

all: $(TARGETS)
clean:
	rm -f $(TARGETS)

vmcw: vmcw.f libs/pde_dp.f he3_funcs.f\
      vmcw_mesh.f vmcw_mon.f vmcw_state.f vmcw_pdecol.f\
      vmcw_cmd.f vmcw_func.f
he3: he3.f he3_funcs.f
vft: vft.f

CFLAGS=-g
FFLAGS=-g

he3dat: he3dat.c he3_funcs.o he3dat.h he3_const.fh
	$(CC) he3dat.c he3_funcs.o -g -lgfortran -lc -lm -o he3dat
#he3_funcs.o: he3_funcs.f
