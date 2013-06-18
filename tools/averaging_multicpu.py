#!/usr/bin/env python2.7

debug = False

import argparse

parser = argparse.ArgumentParser(description='Script to make plots')
parser.add_argument('-s', '--setup', help='Set up jobs', action='store_true')
parser.add_argument('-r', '--run', help='Run (submit) jobs', action='store_true')
parser.add_argument('-a', '--average', help='Do averaging', action='store_true')

args = parser.parse_args()

import sys
import os
import re
import numpy as np

def dist(dist_fn, Ef, BE, BT, T, Nk):
    '''
    takes a distribution function and input parameters.
    returns an array of weights, normalized by default to 1.0
    '''
    arr = dist_fn(Ef, BE, BT, T, Nk)
    #print("arr is "+str(arr))
    #arr /= np.linalg.norm(arr)
    arr /= sum(arr)
    #print("arr is "+str(arr))
    return arr


def fdd(Ef, BE, BT, T, Nk):
    '''
    returns a Fermi-Dirac distribution
    '''
    arr = np.zeros(Nk)
    for i in range(Nk):
        E = (BT - BE)/(Nk-1)*i + BE
        arr[i] = 1/(1 + np.exp((E - Ef)*3.185e5/T))
    return arr


def read_param(paramname, filename):
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
            param = float(re.split('[= \t\n]',line.strip())[1])
            f.close()
            return param
    f.close()
    print "ERROR [read_param]: didn't find %s" % paramname
    return 0


def change_param(paramname, filename, paramvalue):
    '''
    changes a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    output = ''
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip()[0:len(paramname)] == paramname:
                line = re.sub('(?<==)[0-9\.]*', paramvalue, line)
            output += line
    with open(filename, 'w') as f:
        f.write(output)


def averageFile(fileName, w, isTimeDep=True):
    '''
    This function averages the contents of a file, which is expected to be
    in all folders in the directory 'avg/'.
    '''
    # get dimensions of file
    currentFile = 'avg/1/'+fileName
    data = np.loadtxt(currentFile)

    # create output data
    if (len(data.shape) > 1):
        outputData = np.zeros((data.shape[0],data.shape[1]))
    else:
        outputData = np.zeros((data.shape))

    # for each starting state
    for i in range(len(w)):
        # load data file
        currentFile = 'avg/'+str(i+1)+'/'+fileName
        data = np.loadtxt(currentFile)
        # increment weighted average of data
        outputData += w[i]*data

    # if time-dep, copy times rather than keeping average
    if (isTimeDep):
        outputData[:,0] = data[:,0]

    if (debug):
        print 'output data is ', outputData
        print 'output data shape is ', outputData.shape
    # output to file
    if (np.rank(outputData) == 0):
        with open('avg/avg_%s' % fileName, 'w') as f:
            f.write('%f\n' % float(outputData))
    else:
        np.savetxt('avg/avg_%s' % fileName, outputData, '%-.7g')


def do_setup(w):
    '''
    This function sets up input files for each job to be run
    '''
    # for each starting state
    for i in range(len(w)):
        runDir = 'avg/'+str(i+1)
        # make new directory
        os.system('mkdir -p %s' % runDir)
        # copy inputs
        os.system('cp -r ins/ %s' % runDir)
        # change inputs
        change_param('Nk_first', runDir+'/ins/parameters.in', str(i+1))
        change_param('Nk_final', runDir+'/ins/parameters.in', str(i+1))


def do_runs(w):
    '''
    This function submits the jobs created by do_setup
    '''
    # for each starting state
    for i in range(len(w)):
        runDir = 'avg/'+str(i+1)
        # cd to the directory
        os.chdir(runDir)
        # run hoagie with command "$maindir/bin/dynamix && rm -rf ins/"
        os.system('hoagie -k -l -e $(pwd)/../../bin/dynamix')
        # go back to the parent directory
        os.chdir('../..')


def do_averaging(w):
    '''
    This function averages the data in certain output files
    '''
    averageFile('tcprob.out', w)
    averageFile('ms_est.out', w)
    averageFile('ms_est_tot.out', w)
    averageFile('longtime_from_singlestate.out', w, False)
    averageFile('longtime_from_singlestate_sum.out', w, False)


## parameters of FDD distribution
infile = 'ins/parameters.in'
Ef = -0.01				# Fermi level in Hartree
BE = read_param('k_bandedge', infile)	# band edge in Hartree
BT = read_param('k_bandtop', infile)	# band top in Hartree
T = read_param('temp', infile)		# temperature in Kelvin
Nk = int(read_param('Nk', infile))	# number of states in bulk
timesteps = int(read_param('numOutputSteps', infile)) + 1 # timesteps

# array of weights for each run
w = dist(fdd, Ef, BE, BT, T, Nk)
if (debug):
    print 'sum of distribution w is ', sum(w)

if (args.setup):
    print 'doing setup...'
    do_setup(w)
    print 'done with setup!'

if (args.run):
    print 'submitting jobs...'
    do_runs(w)
    print 'done submitting jobs!'

if (args.average):
    print 'averaging data...'
    do_averaging(w)
    print 'done averaging data!'
