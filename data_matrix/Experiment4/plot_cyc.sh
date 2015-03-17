#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'
set output 'cyc.png'
set title "Cyclecount"

set xlabel "N"
set ylabel "Cyclecount"

set key left top

plot "naive_cyc.dat" title 'naive', \
     "naive2_cyc.dat" title 'naive2', \
     "naive2_par2_cyc.dat" title 'naive2-par2', \
     "naive2_par4_cyc.dat" title 'naive2-par4', \
     "col_cyc.dat" title 'A-trans', \
     "row_cyc.dat" linetype 7 title 'B-trans'
