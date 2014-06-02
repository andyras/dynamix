#!/usr/bin/env python2.7

'''
This script processes data files output from dynamix, with an eye for
observables relevant for the torsion+ET problem.
'''

import argparse

parser = argparse.ArgumentParser(description='get the input file')
parser.add_argument('inputFile', help='The input file to this script, output from dynamix')

args = parser.parse_args()

import numpy as np
import os

def writeOutput(fileName, outputArray, fmt='%12.10e'):
    '''
    This function checks for the outs/ directory and outputs there, otherwise
    it dumps to ./
    '''
    if (os.path.isdir('outs/')):
        outPath = 'outs/'+fileName
    else:
        outPath = fileName

    # have to treat scalars differently
    if (outputArray.ndim == 0):
        with open(outPath, 'w') as o:
            fmt += '\n'
            o.write(fmt % outputArray)
    else:
        np.savetxt(fileName, outputArray, fmt=fmt)


input = np.loadtxt(args.inputFile)

dt = input[1,0] - input[0,0]

# compute integral of population
popInt = np.array(np.trapz(input[:,0], x=input[:,1]))
writeOutput('torsion_popInt.out', popInt)

# compute integral of population after peak

# integral while population is above a cutoff

# time above cutoff

# just the peak value