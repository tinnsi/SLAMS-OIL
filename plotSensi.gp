set terminal postscript eps enhanced "Helvetica" 14 lw 2
set output 'figureOil2.ps'
set size 1,0.5
set multiplot
set key bottom right

set size 0.5,0.5
set key bottom right
set origin 0,0
#set logscale x
set xtics 500
set label 'a' at 100,130
#set ylabel "non-petro-carbon Flux to Seafloor [mgm^{-2}day^{-1}]"
set ylabel "OrgC Flux to Seafloor [mgm^{-2}day^{-1}]"
set xlabel 'Primary Production [mgCm^{-2}day^{-1}]'#'{/Symbol a}_{oil}'
plot 'norm3' u ($1*1e3):($3*1e3) t '3 sizes' w lp, 'disp3' u ($1*1e3):($3*1e3) t '1 size' w lp

unset label
set origin 0.5,0
set nokey
set ylabel "Oil Flux to Seafloor [mgm^{-2}day^{-1}]"
unset logscale y
set yrange [0:]
set label 'b' at 100,220
plot 'norm3' u ($1*1e3):($4*1e3) w lp, 'disp3' u ($1*1e3):($4*1e3) w lp


