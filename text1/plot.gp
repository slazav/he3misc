#!/usr/bin/gnuplot

splot "data.txt" using 1:2:3 with lines
pause -1
splot "data.txt" using 1:2:4 with lines
pause -1
