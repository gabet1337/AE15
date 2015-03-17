#!/usr/bin/gnuplot
set terminal pngcairo enhanced font 'Verdana,12'
set output 'l2_tca.png'
set title "L2 Cache Accesses"

set xlabel "N"
set ylabel "Cache Accesses"

set key left top

plot "naive_l2.dat" title 'naive', \
     "naive2_l2.dat" title 'naive2', \
     "naive2_par2_l2.dat" title 'naive2-par2', \
     "naive2_par4_l2.dat" title 'naive2-par4', \
     "col_l2.dat" title 'A-trans', \
     "row_l2.dat" linetype 7 title 'B-trans'
