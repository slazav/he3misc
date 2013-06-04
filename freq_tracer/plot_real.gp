#!/usr/bin/gnuplot

set terminal x11
set style data points

#set terminal png
#set output "img_amp.png"

#  "data_tmp.dat" using 1:2 with lines pt 7,\
#  "data_flt.dat" using 1:2 with lines title "filter" lc 2 pt 7,\

plot [0:5] \
  "data_tmp.dat" using 1:2 with lines pt 7,\
  "data_ph.dat" using 1:2 with lines title "amp" lc 6 pt 7,\
  "data_res1.dat" using 1:3 title "calculated amplitude" lc 3 pt 7,\
  "data_res1.dat" using 1:4 title "error" lc 5 pt 7,\
  "data_fin.dat"  using 1:3 title "freq" lc 4,\
  0 title "" lc 0
pause -1

#set output "img_fre.png"

plot [][0:2200] \
  "data_res1.dat" using 1:2 title "calculated freq" lc 3 pt 7,\
  "20130517_47074.txt" using 1:3 pt 7
pause -1

plot [0:2200] \
  "data_res1.dat" using 2:3 with points lc 3 pt 7,\
  "20130517_47074.txt" using 3:2
pause -1
