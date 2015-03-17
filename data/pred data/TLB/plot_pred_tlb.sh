#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,15'
set output 'test.png'
set title "TLB data misses"
set xlabel "input size"
set ylabel "Misses"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "pred_tlb_lb.dat" title 'lower bound', \
     "pred_tlb_bs.dat" title 'binary search', \
     "pred_tlb_dfs.dat" title 'dfs', \
     "pred_tlb_bfs.dat" title 'bfs', \
     "pred_tlb_inorder.dat" title 'inorder'
