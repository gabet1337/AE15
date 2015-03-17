#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'

set output 'Runtime_1000.png'
set title "Running Time N = 1024"
set xlabel "n"
set ylabel "Runtime (sec.)"
set key left top

plot "rec_runtime_cutoff.dat" title 'cutoff-n'

set terminal pngcairo enhanced font 'Verdana,12'

set output 'Runtime_1000_l2.png'
set title "L2 TCM, N = 1024"
set xlabel "n"
set ylabel "Total Cache Misses"
set key left top

plot "rec_l2_tcm_cutoff.dat" title 'cutoff-n'

set terminal pngcairo enhanced font 'Verdana,12'

set output 'Runtime_1000_tlb.png'
set title "TLB, N = 1024"
set xlabel "n"
set ylabel "Total Cache Misses"
set key left top

plot "rec_tbl_cutoff.dat" title 'cutoff-n'

set terminal pngcairo enhanced font 'Verdana,12'

set output 'Runtime_1000_filter.png'
set title "Running Time N = 1024"
set xlabel "n"
set ylabel "Runtime (sec.)"
set key left top

plot "good_rec_runtime_cutoff_filter.dat" title 'cutoff-n'
