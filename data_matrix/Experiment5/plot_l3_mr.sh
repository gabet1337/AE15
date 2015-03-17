#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'
set output 'cyc.png'
set title "L2 Missratio"

set xlabel "N"
set ylabel "L2 Missratio"

set key left top

plot "naive_l3_mr.dat" title 'naive', \
     "naive2_l3_mr.dat" title 'naive2', \
     "naive_par2_mr.dat" title 'naive2-par2', \
     "naive_par4_mr.dat" title 'naive2-par4', \
     "col_mr.dat" title 'A-trans', \
     "row_mr.dat" linetype 7 title 'B-trans'
