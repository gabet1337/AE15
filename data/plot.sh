#!/usr/bin/gnuplot

set terminal postscript eps enhanced color "Helvetica" 20
set output 'test.png'

set xlabel "input size"
set ylabel "miss ratio (MSP/CN)"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "INO.dat" title 'inorder', \
     "BFS.dat" title 'bfs', \
     "DFS.dat" title 'dfs'
