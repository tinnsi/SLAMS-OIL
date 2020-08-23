set terminal postscript eps "Helvetica" 16 lw 3
set output 'figureCumuOil.ps'
set multiplot
set key top left
set size 1,1

set size 0.48,0.48

set yrange [0:80]
#set logscale x
set origin 0,0.5
set xlabel 'Radius [mm]'
#set ylabel 'orgC flux [%PP]'
#set title '2x, PP=200'
set ylabel 'orgC flux [mg/m^2day]'
#units of Size, Dens, Velo, Oil is in mol/m2/200days og g/m2/200 days
# so we multiply size, dens and velo by 12 to get g/m2/200d and by 5 (1000/200) to get mg/m2/d 
# we multiply by 60 to convert mol/m2/100d to mg/m2/d
# oil is in g/200d so we multiply by 5 to get mg/m2/d.
set label 'a' at 18,10
plot 'Size13674' u ($2*1e-3):($1*60) t '160' w l, 'Size13675' u ($2*1e-3):($1*60) t '320'  w l, 'Size13676' u ($2*1e-3):($1*60) t '640'  w l, 'Size13677' u ($2*1e-3):($1*60) t '1280'  w l, 'Size13678' u ($2*1e-3):($1*60) t '2560'  w l, 'Size13679' u ($2*1e-3):($1*60) t '5120'  w l

set origin 0.5,0.48
set size 0.48,0.5
set nokey
unset label
unset logscale x
set xrange [1.05:1.55]
set xlabel 'Density [g/cm^3]'
#set xrange [1:1.6]
set label 'b' at 1.5,8
plot 'Dens13674' u 2:($1*60) w l, 'Dens13675' u 2:($1*60) w l, 'Dens13676' u 2:($1*60) w l, 'Dens13677' u 2:($1*60) w l, 'Dens13678' u 2:($1*60) w l, 'Dens13679' u 2:($1*60) w l

set origin 0,0
set size 0.48,0.48
unset xrange
unset label 
set label 'c' at 110,10
set xlabel 'Velocity [m/day]'
plot 'Velo13674' u 2:($1*60) w l, 'Velo13675' u 2:($1*60) w l, 'Velo13676' u 2:($1*60) w l, 'Velo13677' u 2:($1*60) w l, 'Velo13678' u 2:($1*60) w l, 'Velo13679' u 2:($1*60) w l

set origin 0.5,0
set xlabel 'Radius [mm]'
set xrange [0:20]
set yrange [0:1000]
unset label 
set logscale y
set ylabel 'oil flux [mg/m^2day]'
#set ylabel 'oil flux [mg/m^2day]'
set label 'd' at 18,0.5
plot 'Oil13674' u ($2*1e-3):($1*5) t '80' w l, 'Oil13675' u ($2*1e-3):($1*5) t '1280' w l, 'Oil13676' u ($2*1e-3):($1*5) w l, 'Oil13677' u ($2*1e-3):($1*5) w l, 'Oil13678' u ($2*1e-3):($1*5) w l, 'Oil13679' u ($2*1e-3):($1*5) w l

