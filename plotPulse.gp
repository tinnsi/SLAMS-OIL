set terminal postscript eps enhanced "Helvetica" 12 #lw 2
set output 'figurePulse.ps'
set size 1,0.5
set multiplot

set style line 1 lt 1 lw 2
set style line 2 lt 4 lw 2
set origin 0,0
set size 0.5,0.5
set logscale y
#set xtics 0.2
set xrange [600:1100]
set title "Oil Flux to Sea Surface"
set ylabel "oil [mg/m^2day]"
set xlabel "time [days]"
set label 'a' at 620,5000
plot 'OutputSLAMS1101/seafloor.dat' u 16:($13*1e3) t 'DwH' w l ls 2,  'OutputSLAMS1105/seafloor.dat' u 16:($13*1e3) t '0.1xDwH' w l ls 1


set origin 0.5,0
unset ylabel
unset label
set title "Oil Flux at Seafloor"
set yrange [0.1:]
set label 'b' at 620,4000
set arrow from 893,0.1 to 893,1
plot 'OutputSLAMS1101/seafloor.dat' u  16:($6*1e3) t 'DwH' w l ls 2,'OutputSLAMS1105/seafloor.dat' u 16:($6*1e3) t '0.1xDwH' w l ls 1, 'stout_HF' u ($1+800):($2/68.8*1000) t 'data' w lp lw 2



