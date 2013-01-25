#!/usr/bin/env python2.7

import re
from math import pi

debug = False

def readParam(paramName, fileName):
    '''reads a parameter (float value) form a file in the format
    paramName=param
    Note: reads the last instance of paramName, since that is what
    will be read in to dynamix'''
    param = None
    with open(fileName, 'r') as f:
        for line in f.readlines():
            if line[0:(len(paramName)+1)] == (paramName + '='):
                param = float(re.split('[= \t\n]', line.strip())[1])
    if (param == None):
        print 'ERROR [readParam]: did not find %s' % paramName
        return 0
    else:
        if (debug):
            print 'parameter %s is %f' % (paramName, param)
        return param

inputFile = '../ins/parameters.in'

bandEdge = readParam('k_bandedge', inputFile)
bandTop = readParam('k_bandtop', inputFile)
N = readParam('Nk', inputFile)

delta = float(bandTop - bandEdge)/(N-1)

tr = 2*pi/delta

print 'Recurrence time due to bulk spacing is at %f time units.' % tr
