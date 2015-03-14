#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'test.png'

set xlabel "input size"
set ylabel "running time (ms)"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "pred_running_time_lb.dat" title 'lower bound', \
     "pred_running_time_bs.dat" title 'binary search', \
     "pred_running_time_dfs.dat" title 'dfs', \
     "pred_running_time_bfs.dat" title 'bfs', \
     "pred_running_time_inorder.dat" title 'inorder'
