#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime_par_naive.png'

set key left top

plot "naive_par_runtime.dat" title 'par', \
     "naive_std_runtime.dat" title 'naive'
