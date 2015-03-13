#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime.png'

set key left top

plot "naive_runtime.dat" title 'naive', \
     "col_runtime.dat" title 'col', \
     "row_runtime.dat" title 'row'
