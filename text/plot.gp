#!/usr/bin/gnuplot

set terminal x11
set style data lines
plot [][0:90]\
  "alpha.dat",\
  "beta.dat"
pause -1