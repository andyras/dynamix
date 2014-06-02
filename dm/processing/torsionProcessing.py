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

dt = input[1,0]
cutoff = 0.02

# compute integral of population
popInt = np.trapz(input[:,1], x=input[:,0])
writeOutput('torsion_popInt.out', popInt)

# compute integral of population after peak
popMaxInd = np.argmax(input[:,1])
popMaxInt = np.trapz(input[popMaxInd:,1], x=input[popMaxInd:,0])
writeOutput('torsion_popMaxInt.out', popMaxInt)

# integral while population is above a cutoff
aboveInd = input[:,1] > cutoff
aboveInt = np.trapz(input[aboveInd,1], x=input[aboveInd,0])
writeOutput('torsion_aboveInt.out', aboveInt)

# time above cutoff
aboveIntTime = sum(aboveInd)*dt
writeOutput('torsion_aboveIntTime', aboveIntTime)

# just the peak value
popMax = np.max(input[:,1])
writeOutput('torsion_popMax.out', popMax)