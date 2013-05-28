#!/usr/bin/gnuplot

set terminal x11
set style data linespoints

#set terminal png
#set output "img_amp.png"

plot [0:0.0011] \
  "data_test.dat" using 1:2 title "signal" lc 1,\
  "data_flt.dat" using 1:2 title "filter -- coordinate"lc 2 pt 7,\
  "data_flt.dat" using 1:3 title "filter -- velocity" lc 3 pt 7,\
  "data_flt.dat" using 1:4 title "damping" lc 4 pt 7,\
  0 title "" lc 0
pause -1

