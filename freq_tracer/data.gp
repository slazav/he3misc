#!/usr/bin/gnuplot

set terminal x11
set style data lines

#set terminal png
#set output "img_amp.png"

plot [] \
  "data1.dat" using 1:2 title "signal" lc 1,\
  "data2.dat" using 1:2 with linespoints lc 2 pt 7,\
  "data3.dat" using 1:4 with points title "calculated amplitude" lc 3 pt 7,\
  "data1.dat" using 1:4 title "amplitude" lc 0,\
  "data1.dat" using 1:5 title "noise" lc 0,\
  "data3.dat" using 1:5 title "error" lc 2,\
  0 title "" lc 0
pause -1

#set output "img_fre.png"

plot [][34500:36000] \
  "data1.dat" using 1:3 title "frequency" lc 1,\
  "data3.dat" using 1:3 with linespoints title "calculated freq" lc 3 pt 7
pause -1
