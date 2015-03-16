#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,13'
set output 'radix_single_TLB.png'
set title "TLB misses single core"
set xlabel "input size"
set ylabel "Miss count"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "radix2_single_TLB_DM.dat" title 'radix 2', \
    "radix4_single_TLB_DM.dat" title 'radix 4', \
    "radix8_single_TLB_DM.dat" title 'radix 8', \
    "radix16_single_TLB_DM.dat" title 'radix 16'

set output 'radix_multi_TLB.png'
set title "TLB misses multi core"
set ylabel "Miss count"

plot "radix2_multi_TLB_DM.dat" title 'radix 2', \
    "radix4_multi_TLB_DM.dat" title 'radix 4', \
    "radix8_multi_TLB_DM.dat" title 'radix 8', \
    "radix16_multi_TLB_DM.dat" title 'radix 16'

