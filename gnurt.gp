
reset

set size 1,1
set origin 0,0
set logscale y
set key at .75,6

set xrange [.6:1.2]
set yrange [.01:10]

set multiplot
set size 1,1./3.
set origin 0,2./3.
set label 1 at 1,.1 '1D' font ",20"

plot 'off_rt.dat' u 1:3 w l t 'Density', '' u 1:6 w l t '% Ejecta'

unset key
set origin 0,1./3.
set label 1 '2D'

plot '2d_rt.dat' w l, '' u 1:10 w l

set origin 0,0
set label 1 '1D + Model'

plot 'on_rt.dat' u 1:3 w l, '' u 1:6 w l


unset multiplot
