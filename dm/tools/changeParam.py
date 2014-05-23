#!/usr/bin/env python2.7
'''
This script changes parameters in a file of the format
---------------------------------------------
paramname=paramvalue # comments and such here
---------------------------------------------
'''

import re
import argparse

debug = False

parser = argparse.ArgumentParser(description='script to change input parameters')
parser.add_argument('paramname', help='Parameter to be changed', type=str)
parser.add_argument('paramvalue', help='Parameter\'s new value', type=str)
parser.add_argument('-filename', '-f', help='File name', type=str, default='ins/parameters.in')

args = parser.parse_args()

if (debug):
    print 'Arguments are:', args

def change_param(paramname, paramvalue, filename):
    '''
    changes a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    foundParam = False
    output = ''
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if (line.strip()[0:len(paramname)+1] == paramname + '='):
                foundParam = True
                line = re.sub('(?<==)[0-9\.\+\-e]*', paramvalue, line)
            output += line
    if (foundParam):
        with open(filename, 'w') as f:
            f.write(output)
    else:
        print 'Warning [%s]: parameter %s not found' % (__file__, paramname)

change_param(args.paramname, args.paramvalue, args.filename)
