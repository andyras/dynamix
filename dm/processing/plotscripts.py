import os

def pdos_plot():
    plotfile = 'pdos.plt'
    with open(plotfile, 'w') as f:
        f.write("#!/usr/bin/env gnuplot\n")
        f.write("set terminal postscript enhanced color size 20cm,10cm font 'Courier-Bold,12'\n")
        f.write("set output 'pdos.eps'\n")
        f.write("set xlabel 'Energy (E_h)'\n")
        f.write("set key outside right\n")
        f.write("set tics scale 0\n")
        f.write("set yr [0:1]\n")
        f.write("set ytics ('0' 0, '1' 1)\n")
        f.write("plot 'bu_proj.out' w im lw 2 t 'Bulk', \\\n")
        f.write("     'qd_proj.out' w im lw 2 lt 1 lc rgb 'green' t 'QD', \\\n")
        f.write("     'br_proj.out' w im lw 2 lt 1 lc rgb 'blue' t 'Bridge'\n")
    os.system('chmod 755 '+plotfile)
    os.system('./'+plotfile)


def pdos_plot_stack():
    plotfile = 'pdos_stack.plt'
    with open(plotfile, 'w') as f:
        f.write("#!/usr/bin/env gnuplot\n")
        f.write("set terminal postscript enhanced color size 20cm,14cm font 'Courier-Bold,12'\n")
        f.write("set output 'pdos_stack.eps'\n")
        f.write("set key outside right\n")
        f.write("set tics scale 0\n")
        f.write("set yr [0:1]\n")
        f.write("set ytics ('0' 0, '1' 1)\n")
        f.write("set rmargin 18\n")
        f.write("set multiplot layout 4,1\n")
        f.write("plot 'bu_proj.out' w im lw 2 lc rgb 'red' t 'Bulk'\n")
        f.write("plot 'qd_proj.out' w im lw 2 lc rgb 'green' t 'QD'\n")
        f.write("plot 'br_proj.out' w im lw 2 lc rgb 'blue' t 'Bridge'\n")
        f.write("set xlabel 'Energy (E_h)'\n")
        f.write("unset yrange\n")
        f.write("plot 'psi_diag.out' w im lw 2 lc rgb 'black' t 'Psi'\n")
        f.write("unset multiplot\n")
    os.system('chmod 755 '+plotfile)
    os.system('./'+plotfile)
