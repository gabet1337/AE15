#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'
set output 'test.png'
set title "Branch Misprediction ratio"
set xlabel "input size"
set ylabel "Miss Ratio (misses/total)"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key right bottom

plot "pred_br_lb.dat" title 'lower bound', \
     "pred_br_bs.dat" title 'binary search', \
     "pred_br_dfs.dat" title 'dfs', \
     "pred_br_bfs.dat" title 'bfs', \
     "pred_br_inorder.dat" title 'inorder'
