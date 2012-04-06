set terminal postscript eps enhanced colour
set output "kprob.eps"
set title "Electron probability density in bulk conduction band"
set border 0
set ytics scale 0
set xtics scale 0
set ylabel "Energy above band edge \(a.u.\)"
set xlabel "Time (a.u.)"
set xrange [0:20000]
set yrange [0.00:.00735004225]
unset key
unset colorbox
set pm3d map
set palette model XYZ functions gray**0.45, gray**1.0, gray**1.0
splot './outs/kprobs_gnuplot.out'
