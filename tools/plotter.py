#!/usr/bin/env python2.7

import os
import sys
import glob

def gnuplot_all(path='../figures'):
    '''
    This function executes all the plot files in the specified directory
    '''
    if (path != '../figures'):
        if not (os.path.isdir(path)):
            path='./figures'
    if not (os.path.isdir(path)):
        print 'BOGUS! Path not found.'
        sys.exit()

    # go to directory and execute each plot file
    os.chdir(path)
    for plot in glob.glob('*.plt'):
        os.system('./'+plot)
