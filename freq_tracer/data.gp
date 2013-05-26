#!/usr/bin/gnuplot

set terminal x11
set style data lines

#set terminal png
#set output "img_amp.png"

plot [] \
  "data_test.dat" using 1:2 title "signal" lc 1,\
  "data_res.dat" using 1:3 with points title "calculated amplitude" lc 3 pt 7,\
  "data_test.dat" using 1:4 title "amplitude" lc 0,\
  "data_test.dat" using 1:5 title "noise" lc 0,\
  "data_res.dat" using 1:4 title "error" lc 2,\
  "data_res.dat" using 1:5 title "base" lc 4,\
  0 title "" lc 0
pause -1

#set output "img_fre.png"

plot [][34500:36000] \
  "data_test.dat" using 1:3 title "frequency" lc 1,\
  "data_res.dat" using 1:2 with linespoints title "calculated freq" lc 3 pt 7
pause -1
