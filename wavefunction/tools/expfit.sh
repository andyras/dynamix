#!/bin/bash

mkdir -p figures

gnuplot 2> fit.log << EOF
set terminal postscript enhanced color lw 2 size 6,4 font 'Arial-Bold,12'
set output '/dev/null'
plot '$1'
fitfn(x) = exp(-k*x)
fit fitfn(x) '$1' using (\$1/GPVAL_X_MAX):2 via k
k = k/GPVAL_X_MAX
k_au = k
k_si = k_au*6.57968*(10**15)
t_au = 1/k_au
t_si = 1/k_si
set output 'figures/expfit.eps'
set key rmargin
set rmargin at screen 0.65
set ylabel 'Population (a.u.)'
set xlabel 'Time (a.u.)'
set tics nomirror scale 0.5
set label "Rate is\n%.2e",k_au," (a.u.)\n%.2e",k_si," (Hz)\n\n1/e time is\n%.2e",t_au," (a.u.)\n%.2e",t_si," (s)" at graph 1.35,0.6 center
plot '$1' t 'Population', fitfn(x) t 'Exp. Decay Fit'
EOF

#rm -f fit.log
