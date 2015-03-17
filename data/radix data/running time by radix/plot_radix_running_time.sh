#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,15'
set output 'radix_single_running_time.png'
set title "Running time single core"
set xlabel "input size"
set ylabel "running time (ms)"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "radix2_single_running_time.dat" title 'radix 2', \
    "radix4_single_running_time.dat" title 'radix 4', \
    "radix8_single_running_time.dat" title 'radix 8', \
    "radix16_single_running_time.dat" title 'radix 16'

set output 'radix_multi_running_time.png'
set title "Running time multi core"

plot "radix2_multi_running_time.dat" title 'radix 2', \
    "radix4_multi_running_time.dat" title 'radix 4', \
    "radix8_multi_running_time.dat" title 'radix 8', \
    "radix16_multi_running_time.dat" title 'radix 16'

