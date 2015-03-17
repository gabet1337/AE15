#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'l3_rec_tca.png'

set key left top

plot "rec_l3.dat" title 'rec', \
     "row2_l3.dat" title 'row'
