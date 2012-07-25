# Makefile

FC=gfortran
LD=ld

texture: texture.f90 var2d.f90 func_text.f90
	gfortran -o $@ $+

test1: texture.f90 var2d.f90 func_test.f90
	gfortran -g -o $@ $+
