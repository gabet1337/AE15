#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,15'
set output 'test.png'
set title "Number of cycles"
set xlabel "input size"
set ylabel "cycle count"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "pred_tot_cyc_lb.dat" title 'lower bound', \
     "pred_tot_cyc_bs.dat" title 'binary search', \
     "pred_tot_cyc_dfs.dat" title 'dfs', \
     "pred_tot_cyc_bfs.dat" title 'bfs', \
     "pred_tot_cyc_inorder.dat" title 'inorder'
