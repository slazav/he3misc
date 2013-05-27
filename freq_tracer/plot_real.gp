#!/usr/bin/gnuplot

set terminal x11
set style data points

#set terminal png
#set output "img_amp.png"

#  "20130517_47074.osc" using ($0 * 2.5E-5 - 10):((($1>=0)?$1-127:$1+127)*6.25000E-02),\

plot [0:45] \
  "20130517_47074.txt" using 1:2 pt 7,\
  "data_res1.dat" using 1:3 title "calculated amplitude" lc 3 pt 7,\
  0 title "" lc 0
pause -1

#set output "img_fre.png"

plot [0:45][650:] \
  "data_res1.dat" using 1:2 title "calculated freq" lc 3 pt 7,\
  "20130517_47074.txt" using 1:3 pt 7
pause -1

plot [650:] \
  "data_res1.dat" using 2:3 with points lc 3 pt 7,\
  "20130517_47074.txt" using 3:2
pause -1
