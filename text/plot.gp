#!/usr/bin/gnuplot

set terminal x11
set style data lines
plot []\
  "texture.dat" using 1:2 title "alpha",\
  "texture.dat" using 1:3 title "beta",\
  "textureoo.dat" using 1:2 title "alpha",\
  "textureoo.dat" using 1:3 title "beta",\
0
pause -1
