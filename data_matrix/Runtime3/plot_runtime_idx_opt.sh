#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime_idx_opt.png'

set key left top

plot "row_idx_runtime.dat" title 'idx', \
     "row_std_runtime.dat" title 'row'
