#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime_par.png'

set key left top

plot "row_par_runtime.dat" title 'par', \
     "row_std_runtime2.dat" title 'row'
