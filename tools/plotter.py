#!/usr/bin/env python2.7

debug = True

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

def makePlotFile(plotFile):
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
        exec plotFn+"(f)"

def plot_test(f):
    print 'whoo'
    f.write("\n")

def plot_vibprob(f):
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
    f.write("set xlabel 'Time (a.u.)'\n")
    f.write("set ylabel 'Vibrational Population'\n")
    f.write("set title 'Vibrational Populations'\n")
    ###
    f.write("set xr [0:%s]\n" % tout)
    f.write("\n")
    f.write("set lmargin 18\n")
    f.write("set key out left\n")
    f.write("\n")
    f.write("plot '../outs/vibprob.out' t '0'")
    for ii in range(1,n):
        f.write(", \\\n")
        f.write("'' u 1:%d t '%d'" % (ii+2, ii))
    f.write("\n")
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

plots = set_plots()

print plots
for plot in plots:
    makePlotFile(plot)
    print plot
    print getPlotFn(plot)

gnuplot_all()
