set terminal postscript "Helvetica" 13 lw 2
set output 'figureSFSed.ps'

set multiplot
set size 1,1

set size 0.49,0.5
#set xtics 500
#set xtics 30

set xrange [0:2000]
set key top left
set origin 0.01,0
set label 'c' at 1800,2.7
set xlabel 'time [days]'
set ylabel 'Total flux at seafloor [g/m^2day]'
set yrange [0:3]
plot 'OutputSLAM3081/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '10' w l, 'OutputSLAM3082/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '50' w l, 'OutputSLAM3083/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '100' w l, 'OutputSLAMS3084/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '400' w l

set key top left
set origin 0.51,0
unset label
set label 'd' at 1800,2.7
set ylabel 'Total flux at seafloor [g/m^2day]'
plot 'OutputSLAM3085/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '10' w l, 'OutputSLAM3086/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '50' w l, 'OutputSLAM3087/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '100' w l, 'OutputSLAM3088/seafloor.dat' u  16:($1+$2+$3+$4+$5+$6+$7) t '400' w l

set origin 0,0.5
set size 0.5,0.5
set title 'Surface Oil Flux = 1600mg/m^2day'
unset label
set label 'a' at 1800,300
set ylabel 'Oil flux at seafloor [mg/m^2day]'
set logscale y
set yrange [0.01:1000]
unset xlabel
plot 'OutputSLAM3081/seafloor.dat' u  16:($6*1e3) t '10' w l, 'OutputSLAM3082/seafloor.dat' u  16:($6*1e3) t '50' w l, 'OutputSLAM3083/seafloor.dat' u  16:($6*1e3) t '100' w l, 'OutputSLAMS3084/seafloor.dat' u  16:($6*1e3) t '400' w l

set title 'Surface Oil Flux = 400mg/m^2day'
set origin 0.5,0.5
unset label
set label 'b' at 1800,300
set ylabel 'Oil flux at seafloor [mg/m^2day]'
set logscale y
plot 'OutputSLAM3085/seafloor.dat' u  16:($6*1e3) t '10' w l, 'OutputSLAM3086/seafloor.dat' u  16:($6*1e3) t '50' w l, 'OutputSLAM3087/seafloor.dat' u  16:($6*1e3) t '100' w l, 'OutputSLAM3088/seafloor.dat' u  16:($6*1e3) t '400' w l

