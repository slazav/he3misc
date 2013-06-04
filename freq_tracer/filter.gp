#!/usr/bin/gnuplot

set terminal x11
set style data linespoints

#set terminal png
#set output "img_amp.png"

plot [0.0007:0.0009][-3:3] \
  "data_test.dat" using 1:2 title "signal" lc 1,\
  "data_flt.dat" using 1:3 title "filter1" lc 3 pt 7,\
  "data_ph.dat" using 1:2 title "ph" lc 4 pt 7,\
  "data_flt.dat" using 1:4 title "damping" lc 4 pt 7,\
  0 title "" lc 0
pause -1

