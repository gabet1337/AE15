#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime.png'

set xlabel "input size"
set ylabel "runtime (sec)"

set key left top

plot "naive_runtime.dat" title 'naive', \
     "naive2_runtime.dat" title 'naive2', \
     "naive2_par2_runtime.dat" title 'naive2-par2', \
     "naive2_par4_runtime.dat" title 'naive2-par4', \
     "col_runtime.dat" title 'col', \
     "row_runtime.dat" title 'row'
