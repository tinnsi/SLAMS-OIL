set terminal postscript enhanced 
set output 'particle3663.ps'
set size 1,1
set multiplot

set nokey

set origin 0.01,0
set size 0.49,0.5
#set xtics 1000
set xrange [800:1000]
set xlabel 'Time [days]'
set ylabel 'Aggregate Radius [mm]'
set label 'b' at 980,18
plot 'OutputSLAMS3663/particlesFloor.dat' u ($13/2):($2*1e-3)

unset label
set origin 0.0,0.5
set size 0.5,0.5
set ylabel 'oil flux [mg/agg class]'
set label 'a' at 980,160
plot 'OutputSLAMS3663/particlesFloor.dat' u ($13/2):($10*1e3)

unset label
set origin 0.51,0.5
set size 0.49,0.5
set xrange [0:120]
set ylabel 'density [g/cm^3]' offset 0
set xlabel 'velocity [m/day]'
set label 'c' at 110,2.8
set yrange [1:3]
set ytics 1
plot 'OutputSLAMS3663/particlesFloor.dat' u ($3*864):4 t '3663'

unset label
set ytics 2
set origin 0.5,0.0
set size 0.5,0.5
set ylabel 'Aggregate Radius [mm]' offset 1
unset yrange
set label 'd' at 110,18
plot 'OutputSLAMS3663/particlesFloor.dat' u ($3*864):($2*1e-3) t '3663'




