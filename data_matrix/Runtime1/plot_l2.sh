#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'l2_tca.png'

set key left top

plot "naive_l2.dat" title 'naive', \
     "col_l2.dat" title 'col', \
     "row_l2.dat" title 'row'
