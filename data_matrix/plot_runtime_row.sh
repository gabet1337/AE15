#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'runtime.png'

set xlabel "input size"
set ylabel "runtime (sec)"

set key left top

plot "row_runtime_2.dat" title 'row'
