#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'l2_rec_tca.png'

set key left top

plot "rec_l2.dat" title 'col', \
     "row2_l2.dat" title 'row'
