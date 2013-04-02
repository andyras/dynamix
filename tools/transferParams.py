#!/usr/bin/env python2.7
'''
This script applies the parameters in an input file to an output file.
Parameters in the input are of the format
---------------------------------------------
paramname=paramvalue
---------------------------------------------
'''

import re
import argparse

debug = False

parser = argparse.ArgumentParser(description='Script to apply parameters from one file to another')
parser.add_argument('--verbose', '-v', help='Verbose output', action='store_true')
parser.add_argument('--quiet', '-q', help='Quiet output', action='store_true')
parser.add_argument('--dryrun', '-d', help='Dry run: nothing changed in output', action='store_true')
parser.add_argument('--force', '-f', help='Force prepend of parameters to output', action='store_true')
parser.add_argument('-output', '-o', help='Output file name', type=str)
parser.add_argument('-input', '-i', help='Input file name', type=str)

args = parser.parse_args()

if (debug):
    print 'Arguments are:', args

def readParam(paramname, filename):
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

def changeParam(paramname, paramvalue, filename):
    '''
    changes a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    foundParam = False
    output = ''
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if (line.find(paramname) == 0):
                foundParam = True
                line = re.sub('(?<==)[0-9\.\+\-e]*', paramvalue, line)
            output += line
    if (foundParam):
        with open(filename, 'w') as f:
            f.write(output)
    else:
        print 'Warning [%s]: parameter %s not found' % (__file__, paramname)

def getParam(line):
    '''
    gets the name and value of a parameter from a file. The name is assumed to
    be everything before the equals sign, and the value everything after.
    '''
    param = ''
    # ignore commented lines
    if (line[0] == '#'):
        return (param, '')
    else:
        eqpos = line.find('=')
        # no parameter if equals is missing, at beginning or end
        if ((eqpos == -1) or (eqpos == 0) or (eqpos == len(line))):
            return (param, '')
        # return everything before equals sign
        if (debug):
            print "line is %s" % line,
            print "position of equals is %d" % eqpos
            print "param is %s" % line[0:eqpos]
            print "param value is %s" % line[eqpos+1:]
        return (line[0:eqpos], re.split('[= \t\n]',line.strip())[1])

inputFile = args.input
outputFile = args.output

# go through input line by line
with open(inputFile, 'r') as infile:
    lines = infile.readlines()

for line in lines:
    (param, paramValue) = getParam(line)
    # if line has a parameter
    if (not (param == '')):
        oldParamValue = readParam(param, outputFile)
        if (debug):
            print "Old parameter value is %s" % oldParamValue
        # no change if values are the same
        if (oldParamValue == paramValue):
            if (args.verbose):
                print "Not changing %s from %s" % (param, paramValue)
        # parameter not found in output
        elif (oldParamValue == ''):
            if (not args.quiet):
                print "Parameter '%s' not found in output file." % param
            # prepend parameter to file if 'force' option is set
            if (args.force):
                if (not args.quiet):
                    print "prepending '%s' to output file." % line.strip()
                with open(outputFile, 'r') as f: data = f.read()
                with open(outputFile, 'w') as f: f.write(line + data)
        else:
            changeParam(param, paramValue, outputFile)
            if (not args.quiet):
                print "Changing %s to %s" % (param, paramValue)

# for each parameter
