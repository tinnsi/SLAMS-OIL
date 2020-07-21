set terminal postscript eps enhanced "Helvetica" 14 lw 2
set output 'figureSensitivity2.ps'
set size 1,0.5
set multiplot
set nokey

set origin 0,0
set size 0.5,0.5
set yrange [0:300]
set ylabel "Oil flux to Seafloor [mgm^{-2}day^{-1}]" offset 1
set xlabel "{/Symbol a}_{oil} "
set label "a" at 0.15,270
plot 'alpha_oil' u 10:(1e3*$4) w lp

unset label
set origin 0.5,0.02
set size 0.5,0.48
set nokey
#set logscale x
set xlabel "TEP production [mgCm^{-2}day^{-1}]"
set label "b" at 12,270
plot 'TEP_prod' u ($11*1e3):(1e3*$4) w p


