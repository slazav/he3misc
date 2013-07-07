FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

TARGETS=vmcw vft he3
all: $(TARGETS)

###################################

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_pmelt he3_tab he3_tc he3_vm\
        he3_meff he3_flegg he3_d_exp\
        he3_gammaf he3_dnde he3_vf he3_pf\
        he3_yosida he3_z0 he3_f0a he3_susept\
        he3_swvel he3_tau_lt

# additional fitting functions used in libhe3
ADDOBJS=E02AEE E02CBE M01AGE P01AAE X02AAE X04AAE

# Legget equations
LEGG_EQ_OBJS=he3b_legg_rot1d

libhe3.a:\
  $(patsubst %,src/%.o,$(LIBOBJS))\
  $(patsubst %,libs/%.o,$(ADDOBJS))\
  $(patsubst %,legg_eq/%.o,$(LEGG_EQ_OBJS))
	ar rs $@ $+

###################################

clean:
	rm -f $(TARGETS) src/*.o libs/*.o *.o *.a

vmcw: vmcw.f\
      vmcw_mesh.f vmcw_mon.f vmcw_state.f vmcw_pdecol.f\
      vmcw_cmd.f vmcw_func.f libs/pde_dp.f libhe3.a
he3: he3.f libhe3.a
vft: vft.f


CFLAGS=-g
FFLAGS=-g

he3dat: he3dat.c he3_funcs.o he3dat.h he3_const.fh
	$(CC) he3dat.c he3_funcs.o -g -lgfortran -lc -lm -o he3dat
#he3_funcs.o: he3_funcs.f
