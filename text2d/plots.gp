#!/usr/bin/gnuplot

set terminal x11
set style data lines
plot\
  "spec.dat" using 1:2 title "new",\
  "speco.dat" using 1:2 title "old",\
0
pause -1
