#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,15'
set output 'radix_running_time.png'
set title "Running time"
set xlabel "input size"
set ylabel "running time (ms)"
set xtics autofreq
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

set logscale x 2
set format x "2^{%L}"

set key left top

plot "radix_single_running_time.dat" title 'single', \
     "radix_multi_running_time.dat" title 'multi'

set output 'radix_l2_ratio.png'
set title "L2 cache miss ratio"
set ylabel "Miss ratio"

plot "radix_single_L2_ratio.dat" title 'single', \
     "radix_multi_L2_ratio.dat" title 'multi'

set output 'radix_l3_ratio.png'
set title "L3 cache miss ratio"
set ylabel "Miss ratio"

plot "radix_single_L3_ratio.dat" title 'single', \
     "radix_multi_L3_ratio.dat" title 'multi'

set output 'radix_tlb_misses.png'
set title "TLB Misses"
set ylabel "Miss count"

plot "radix_single_TLB_DM.dat" title 'single', \
     "radix_multi_TLB_DM.dat" title 'multi'

set output 'radix_br_msp.png'
set title "Branch Misprediction ratio"
set ylabel "Miss prediction ratio"

plot "radix_single_BR_MSP.dat" title 'single', \
     "radix_multi_BR_MSP.dat" title 'multi'

