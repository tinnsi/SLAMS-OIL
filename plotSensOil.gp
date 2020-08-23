set terminal postscript eps enhanced "Helvetica" 14 lw 2
set output 'figureOilF.ps'
set size 1,1
set multiplot

set size 0.5,0.48
set key top right
set origin 0,0.5
set logscale x
set label 'a' at 12,70
set ylabel "OrgC Flux to Seafloor [mgm^{-2}day^{-1}]"
plot 'oil1000' u ($7*1e3):($3*1e3) t '1000' w lp, 'oil500' u ($7*1e3):($3*1e3) t '500' w lp

unset label
set origin 0.5,0.5
set ylabel "Oil Flux to Seafloor [mgm^{-2}day^{-1}]"
unset logscale y
set arrow from 1600,0 to 1600,100
set yrange [0:]
set label 'b' at 12,400
plot 'oil1000' u ($7*1e3):($4*1e3) t '1000' w lp, 'oil500' u ($7*1e3):($4*1e3) t '500' w lp

unset label
unset arrow
set origin -0.01,0
set key bottom right
set size 0.51,0.52
set xlabel "Oil Flux to Sea Surface [mgm^{-2}day^{-1}]"
set ylabel "oil/(orgC+oil) at seafloor in %" offset 1
set label 'c' at 12,90
plot 'oil1000' u ($7*1e3):($4/($4+$3)*1e2) t '1000' w lp, 'oil500' u ($7*1e3):($4/($4+$3)*1e2) t '500' w lp

unset label
set origin 0.51,0
set size 0.49,0.52
set label 'd' at 12,30
set arrow from 1600,0 to 1600,15
set ylabel "% of Surface Oil at Seafloor" offset 1
#plot 'oil1000' u ($7*1e3):($4/$7*100) t '1000' w lp, 'oil500' u ($7*1e3):($4/$7*100) t '500' w lp
plot 'oilR1000' u ($6*1e3):($5*100) t '1000' w l, 'oilR500' u ($6*1e3):($5*100)t '500'  w l

