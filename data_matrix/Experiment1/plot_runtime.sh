#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'
set output 'runtime.png'
set title "Running time"

set xlabel "N"
set ylabel "running time (sec)"

set key left top

plot "naive_runtime.dat" title 'naive', \
     "naive2_runtime.dat" title 'naive2', \
     "naive2_par2_runtime.dat" title 'naive2-par2', \
     "naive2_par4_runtime.dat" title 'naive2-par4', \
     "col_runtime.dat" title 'A-trans', \
     "row_runtime.dat" linetype 7 title 'B-trans'
