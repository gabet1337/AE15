#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'test.png'

plot "Naive.dat" title 'Naive', \
     "Col.dat" title 'Col', \
     "Row.dat" title 'Row'
