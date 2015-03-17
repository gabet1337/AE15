#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'l3_tca.png'

set key left top

plot "naive_l3.dat" title 'naive', \
     "col_l3.dat" title 'col', \
     "row_l3.dat" title 'row'
