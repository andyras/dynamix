set terminal postscript enhanced colour dashed lw 2 size 6,4 font 'Arial-Bold,10'
set output 'time_and_space.eps'

set multiplot

set style data lines
set style increment user
set style line 1 lt 1 lc rgb 'red'
set style line 2 lt 1 lc rgb 'green'
set style line 3 lt 1 lc rgb 'blue'

set lmargin 24
set key lmargin

set xrange [0:20000]

set ytics nomirror scale 0.5
unset xtics

set size 1, 0.33
set origin 0.0,0.66
plot 'outs/tkprob.out' u 1:2 t 'Bulk' lt 1

set size 1, 0.33
set origin 0.0,0.0
set xlabel 'Time (a.u.)'
set xtics nomirror scale 0.5
plot 'outs/tcprob.out' u 1:2 t 'QD' lt 2
unset xtics
unset xlabel

set size 1, 0.33
set origin 0.0,0.33
set ylabel 'Population (a.u.)'
plot 'outs/tbprob.out' u 1:2 t 'Bridge' lt 3
