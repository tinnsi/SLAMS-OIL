set terminal postscript enhanced solid "Helvetica" 11 lw 2
set output 'figureModVal.ps'

set size 1,0.4
set multiplot
set nokey
set xrange [-0.5:3.5]
set xtics 1

set size 0.38,0.4
set origin 0,0
set title "DWH"
set ylabel 'flux [mg/m^2day]'
set label 'a' at -0.2,400
set logscale y
set yrange [1:1000]
plot 'model_D' u 2:xticlabels(1) w p pt 5, 'giering_D_o' w errorbars

set origin 0.35,0
set title "SEEP"
unset ylabel
unset label
set label 'b' at -0.2,400
set size 0.35,0.4
unset ylabel
#set format y ""
plot 'model_S' u 2:xticlabels(1) pt 5, 'giering_S' w errorbars

set origin 0.66,0
set size 0.35,0.4
unset label
set label 'c' at -0.2,400
unset format y
set title "REF"
plot 'model_r' u 2:xticlabels(1) pt 5, 'giering_R' w errorbars
