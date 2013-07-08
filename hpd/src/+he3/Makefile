LIBNAME=libhe3
MEX=/opt/matlab-R2011b/bin/mex

FUNC_I1 = he3_pmelt he3_tc he3_tab
FUNC_I0 = he3_pa he3_pabn he3_tabn he3_amass

mex: he3_tc.mexglx he3_tab.mexglx he3_pmelt.mexglx

%.mexglx: mexfunc.c ../$(LIBNAME).so
	$(MEX) -glnx86 -DFUNC=$*_\
	   mexfunc.c ../$(LIBNAME).so -output $*

######

%.mexa64: %.o mexfunc.c ../$(LIBNAME).so
	$(MEX) -v -glnxa64 -DFUNC=$*_\
	   mexfunc.c ../$(LIBNAME).so -output $*

clean:
	rm -f *.a *.so *.o *.mexglx
