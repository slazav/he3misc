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
  -Wno-unused-parameter -fno-range-check\
  -I../lib

TARGETS=vmcw vft
all: $(TARGETS)

FC=gfortran

clean:
	rm -f $(TARGETS) *.o

vmcw: vmcw.f\
      vmcw_mesh.f vmcw_mon.f vmcw_state.f vmcw_pdecol.f\
      vmcw_cmd.f vmcw_func.f vmcw_pnm.o ../lib/libs/pde_dp.f ../lib/libhe3.a
vft: vft.f

