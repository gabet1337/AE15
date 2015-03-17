#!/usr/bin/gnuplot

set terminal pngcairo enhanced font 'Verdana,12'

set output 'plot/L2_L3_TCA_runtime.png'
set title "Running Time"
set xlabel "N"
set ylabel "Runtime (sec.)"
set key left top

plot "naive_l2_l3_runtime.dat" title 'naive', \
     "naive2_l2_l3_runtime.dat" title 'naive2', \
     "naive2_par2_l2_l3_runtime.dat" title 'naive2-par2', \
     "naive2_par4_l2_l3_runtime.dat" title 'naive2-par4', \
     "a-trans_l2_l3_runtime.dat" title 'A-trans', \
     "b-trans_l2_l3_runtime.dat" linetype 7 title 'B-trans'

set output 'plot/L2_TCA.png'
set title "L2 TCA"
set xlabel "N"
set ylabel "Cache Accesses"
set key left top

plot "naive_l2_l3_l2_tca.dat" title 'naive', \
     "naive2_l2_l3_l2_tca.dat" title 'naive2', \
     "naive2_par2_l2_l3_l2_tca.dat" title 'naive2-par2', \
     "naive2_par4_l2_l3_l2_tca.dat" title 'naive2-par4', \
     "a-trans_l2_l3_l2_tca.dat" title 'A-trans', \
     "b-trans_l2_l3_l2_tca.dat" linetype 7 title 'B-trans'

set output 'plot/L3_TCA.png'
set title "L3 TCA"
set xlabel "N"
set ylabel "Cache Accesses"
set key left top

plot "naive_l2_l3_l3_tca.dat" title 'naive', \
     "naive2_l2_l3_l3_tca.dat" title 'naive2', \
     "naive2_par2_l2_l3_l3_tca.dat" title 'naive2-par2', \
     "naive2_par4_l2_l3_l3_tca.dat" title 'naive2-par4', \
     "a-trans_l2_l3_l3_tca.dat" title 'A-trans', \
     "b-trans_l2_l3_l3_tca.dat" linetype 7 title 'B-trans'

set output 'plot/TLB_DM.png'
set title "TLB Data Misses"
set xlabel "N"
set ylabel "Misses"
set key left top

plot "naive_tlb_cyc_tlb.dat" title 'naive', \
     "naive2_tlb_cyc_tlb.dat" title 'naive2', \
     "naive2_par2_tlb_cyc_tlb.dat" title 'naive2-par2', \
     "naive2_par4_tlb_cyc_tlb.dat" title 'naive2-par4', \
     "a-trans_tlb_cyc_tlb.dat" title 'A-trans', \
     "b-trans_tlb_cyc_tlb.dat" linetype 7 title 'B-trans'

set output 'plot/CYC.png'
set title "Total Cycle Count"
set xlabel "N"
set ylabel "Cycle Count"
set key left top

plot "naive_tlb_cyc_cyc.dat" title 'naive', \
     "naive2_tlb_cyc_cyc.dat" title 'naive2', \
     "naive2_par2_tlb_cyc_cyc.dat" title 'naive2-par2', \
     "naive2_par4_tlb_cyc_cyc.dat" title 'naive2-par4', \
     "a-trans_tlb_cyc_cyc.dat" title 'A-trans', \
     "b-trans_tlb_cyc_cyc.dat" linetype 7 title 'B-trans'

set output 'plot/L2_CMR.png'
set title "L2 Cache Miss Ratio"
set xlabel "N"
set ylabel "Miss Ratio"
set key left top

plot "naive_l2_miss_ratio.dat" title 'naive', \
     "naive2_l2_miss_ratio.dat" title 'naive2', \
     "naive2_l2_miss_ratio.dat" title 'naive2-par2', \
     "naive2_l2_miss_ratio.dat" title 'naive2-par4', \
     "a-trans_l2_miss_ratio.dat" title 'A-trans', \
     "b-trans_l2_miss_ratio.dat" linetype 7 title 'B-trans'

set output 'plot/L3_CMR.png'
set title "L3 Cache Miss Ratio"
set xlabel "N"
set ylabel "Miss Ratio"
set key left top

plot "naive_l3_miss_ratio.dat" title 'naive', \
     "naive2_l3_miss_ratio.dat" title 'naive2', \
     "naive2_l3_miss_ratio.dat" title 'naive2-par2', \
     "naive2_l3_miss_ratio.dat" title 'naive2-par4', \
     "a-trans_l3_miss_ratio.dat" title 'A-trans', \
     "b-trans_l3_miss_ratio.dat" linetype 7 title 'B-trans'
