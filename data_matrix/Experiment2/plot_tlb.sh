#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'
set output 'tlb.png'
set title "TLB datamisses"

set xlabel "N"
set ylabel "TLB datamisses"

set key left top

plot "naive_tlb.dat" title 'naive', \
     "naive2_tlb.dat" title 'naive2', \
     "naive2_par2_tlb.dat" title 'naive2-par2', \
     "naive2_par4_tlb.dat" title 'naive2-par4', \
     "col_tlb.dat" title 'A-trans', \
     "row_tlb.dat" linetype 7 title 'B-trans'
