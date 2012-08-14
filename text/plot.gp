#!/usr/bin/gnuplot

set terminal x11
set style data lines
plot [][0:90]\
  "texture.dat" using 1:2 title "alpha",\
  "texture.dat" using 1:3 title "beta"
pause -1

plot \
  "alpha.dat"
pause -1