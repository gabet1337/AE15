#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime2.png'

set key left top

plot "rec_runtime.dat" title 'rec', \
     "row2_runtime.dat" title 'row'
