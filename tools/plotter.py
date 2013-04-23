#!/usr/bin/env python2.7

debug = True

import argparse

parser = argparse.ArgumentParser(description='Script to make plots')
parser.add_argument('--plot', '-p', help='Make plots with gnuplot (default is only to create scripts)', action='store_true')
parser.add_argument('--au', '-a', help='Atomic time units (default fs)', action='store_true')

args = parser.parse_args()

import os
import re
import sys
import glob

def readParam(paramname, filename='../ins/parameters.in'):
    '''
    reads a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    f = open(filename)
    lines = f.readlines()
    for line in lines:
        if line[0:len(paramname)] == paramname:
            #print "found %s" % paramname
            #print line
            param = re.split('[= \t\n]',line.strip())[1]
            f.close()
            return param
    f.close()
    if (args.verbose):
        print "WARNING [readParam]: didn't find %s" % paramname
    return ''

def set_plots():
    '''
    Set the plots array to contain the names of all the plot files to be made
    '''
    plots = []
    plots.append("test.plt")
    plots.append("vibprob.plt")
    plots.append("vibprob_subsystem.plt")

    return plots

def gnuplot_all():
    '''
    This function executes all the plot files in the specified directory
    '''

    # execute each plot file
    for plot in glob.glob('*.plt'):
        if (debug):
            print("Plotting %s..." % plot)
        os.system('gnuplot ./'+plot)

def getPlotFn(plotFile):
    '''
    Returns the name (as a string) of the function corresponding to a certain
    plot file
    '''

    # get file name and extension
    (fileName, fileExt) = os.path.splitext(plotFile)
    
    if (fileExt not in ['.plt', '.gp', '.plot', '.gnuplot']):
        print("WARNING: %s is not a normal gnuplot file extension." % fileExt)

    plotFn = "plot_"+fileName
    if (debug):
        print plotFn

    return plotFn

def makePlotFile(plotFile, params):
    '''
    Executes the plot function corresponding to a certain plot file
    '''

    plotFn = getPlotFn(plotFile)

    if not (plotFn in globals()):
        print("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
        print("The function %s to plot %s is not defined." % (plotFn, plotFile))
        print globals()
        sys.exit()

    if (debug):
        print("Executing %s." % plotFn)

    with open(plotFile, 'w') as f:
        exec plotFn+"(f,params)"

def plot_test(f, p):
    print 'whoo'
    f.write("\n")

def plot_vibprob(f, p):
    n = int(readParam("N_vib"))
    tout = readParam("tout")

    f.write("\n")
    f.write("#!/usr/bin/env gnuplot\n")
    f.write("\n")
    f.write("reset\n")
    f.write("set terminal pdfcairo dashed enhanced color size 4,3 font 'FreeMono,12' lw 4 dl 1\n")
    f.write("set style data lines\n")
    f.write("set output 'vibprob.pdf'\n")
    f.write("\n")
    if (p["units"] == 'au'):
        f.write("set xlabel 'Time (a.u.)'\n")
    else:
        f.write("set xlabel 'Time (fs)'\n")
    f.write("set ylabel 'Vibrational Population'\n")
    f.write("set title 'Vibrational Populations'\n")
    if (p["units"] == 'au'):
        f.write("set xr [0:%s]\n" % tout)
    else:
        f.write("set xr [0:%s/41.3414]\n" % tout)
    f.write("\n")
    f.write("set lmargin 18\n")
    f.write("set key out left\n")
    f.write("\n")
    if (p["units"] == 'au'):
        f.write("plot '../outs/vibprob.out' t '0'")
    else:
        f.write("plot '../outs/vibprob.out' u ($1/41.3414):2 t '0'")
    for ii in range(1,n):
        f.write(", \\\n")
        if (p["units"] == 'au'):
            f.write("'' u 1:%d t '%d'" % (ii+2, ii))
        else:
            f.write("'' u ($1/41.3414):%d t '%d'" % (ii+2, ii))
    f.write("\n")
    f.write("\n")
    f.write("reset\n")

def plot_vibprob_subsystem(f, p):
    n = int(readParam("N_vib"))
    tout = readParam("tout")

    f.write("#!/usr/bin/env gnuplot\n")
    f.write("\n")
    f.write("reset\n")
    f.write("\n")
    f.write("set terminal pdfcairo dashed enhanced color size 8,5 font 'FreeMono,12' lw 4 dl 1\n")
    f.write("set style data lines\n")
    f.write("set output 'vibprob_subsystem.pdf'\n")
    f.write("\n")
    if (p["units"] == 'au'):
        f.write("set xlabel 'Time (a.u.)'\n")
    else:
        f.write("set xlabel 'Time (fs)'\n")
    if (p["units"] == 'au'):
        f.write("set xr [0:%s]\n" % tout)
    else:
        f.write("set xr [0:%s/41.3414]\n" % tout)
    f.write("\n")
    f.write("set lmargin 18\n")
    f.write("set key out left font ',6'\n")
    f.write("\n")
    f.write("set multiplot layout 3,1 title 'Subsystem Vibrational Populations'\n")
    f.write("\n")
    f.write("set ylabel 'Bulk'\n")
    if (p["units"] == 'au'):
        f.write("plot '../outs/vibprob_bu.out' t '0'")
    else:
        f.write("plot '../outs/vibprob_bu.out' u ($1/41.3414):2 t '0'")
    for ii in range(1,n):
        f.write(", \\\n")
        if (p["units"] == 'au'):
            f.write("'' u 1:%d t '%d'" % (ii+2, ii))
        else:
            f.write("'' u ($1/41.3414):%d t '%d'" % (ii+2, ii))
    f.write("\n\n")
    f.write("set ylabel 'Bridge'\n")
    if (p["units"] == 'au'):
        f.write("plot '../outs/vibprob_br.out' t '0'")
    else:
        f.write("plot '../outs/vibprob_br.out' u ($1/41.3414):2 t '0'")
    for ii in range(1,n):
        f.write(", \\\n")
        if (p["units"] == 'au'):
            f.write("'' u 1:%d t '%d'" % (ii+2, ii))
        else:
            f.write("'' u ($1/41.3414):%d t '%d'" % (ii+2, ii))
    f.write("\n\n")
    f.write("set ylabel 'QD'\n")
    if (p["units"] == 'au'):
        f.write("plot '../outs/vibprob_qd.out' t '0'")
    else:
        f.write("plot '../outs/vibprob_qd.out' u ($1/41.3414):2 t '0'")
    for ii in range(1,n):
        f.write(", \\\n")
        if (p["units"] == 'au'):
            f.write("'' u 1:%d t '%d'" % (ii+2, ii))
        else:
            f.write("'' u ($1/41.3414):%d t '%d'" % (ii+2, ii))
    f.write("\n\n")
    f.write("unset multiplot\n")
    f.write("\n")
    f.write("reset\n")

# This script will do everything in the directory where .plt files and outputs live
if not ('path' in locals()):
    # try a couple of defaults if nothing was specified
    path='../figures'
    if not (os.path.isdir(path)):
        path='./figures'

if not (os.path.isdir(path)):
    print 'BOGUS! Path not found.'
    sys.exit()

# go to the directory since we know it exists.
os.chdir(path)

# get list of plots to make
plots = set_plots()

# set up parameters
params = {}
if (args.au):
    params['units'] = 'au'
else:
    params['units'] = 'fs'

print plots
for plot in plots:
    makePlotFile(plot, params)
    print plot
    print getPlotFn(plot)

gnuplot_all()
