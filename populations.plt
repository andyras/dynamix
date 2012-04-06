set terminal postscript eps enhanced colour dashed lw 2 size 6,4 font 'Arial-Bold,20'
set output 'populations.eps'
set multiplot
set style data lines
set style line 1 lw 2 lt 1 lc rgb 'red'
set style line 2 lw 2 lt 1 lc rgb 'green'
set style line 3 lw 2 lt 1 lc rgb 'blue'
set style line 4 lw 2 lt 1 lc rgb 'purple'
set style line 5 lw 2 lt 1 lc rgb 'dark-goldenrod'
set style increment user
set key lmargin
set lmargin 24
set border 10
set ytics scale 0.5
set ylabel 'Population (a.u.)'
set xrange [0:483.7765]
set xtics nomirror scale 0.5
set size 1, 0.7; set origin 0.0,0.3
plot 'outs/tkprob.out' using ($1/41.34):2 title 'k states' lt 1, \
'outs/tcprob.out' using ($1/41.34):2 title 'c states' lt 2
set size 1, 0.3
set origin 0.0,0.0
set format y '%.2e'
set ytics ('0' 0 ,sprintf('%.1e',0.0009421) 0.0009421)
set xlabel 'Time (fs)'
set ylabel 'Bridge Pop. (a.u.)'
plot 'outs/bprobs.out' using ($1/41.34):2 title 'b1' lt 3