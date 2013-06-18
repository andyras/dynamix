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
    arr /= np.linalg.norm(arr)
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


def do_runs2(infile, w, timesteps, redo=False):
    # create output array variable
    output_avg = np.zeros((timesteps, 2))
    if redo or not os.path.isdir('avg'):
        # make directory for averaging
        os.system('rm -rf avg')
        os.mkdir('avg')
    # for each starting state
    for i in range(len(w)):
        run_file = 'avg/tcprob_'+str(i+1)+'.out'
        # if current weight is > 0.001 of largest weight
        if w[i] >= 0.001*w.max():
            # if redo flag is on and file doesn't exist
            if redo or not os.path.isfile(run_file):
                print '\n\n\nwhooooo\n\n\n'
                # change parameters in ins/parameters.sh
                change_param('Nk_first', infile, str(i+1))
                change_param('Nk_final', infile, str(i+1))
                # run total_dynamix
                # os.system('./total_dynamix')
                os.system('../../dynamix')
                # cp output(s) to folder for averaging
                os.system('cp tcprob.out '+run_file)
            # load values from output
            data = np.loadtxt('avg/tcprob_'+str(i+1)+'.out')
            # add weighted contribution to output
            print('shape of tcprob data file '+str(i)+' is '+str(data.shape))
            for j in range(data.shape[0]):
                output_avg[j,1] += data[j,1]*w[i]
        # add times to output
        for j in range(data.shape[0]):
            output_avg[j,0] = data[j,0]
    # write output to file
    np.savetxt('avg/tcprob_avg.out', output_avg, '%-.7g')


def averageFile(fileName, w):
    # TODO add flag for whether time-dep property
    '''
    This function averages the contents of a file, which is expected to be
    in all folders in the directory 'avg/'.
    '''
    # get dimensions of file
    currentFile = 'avg/1/'+fileName
    data = np.loadtxt(currentFile)
    nr = data.shape[0]  # number of rows
    nc = data.shape[1]  # number of cols

    # create output data
    outputData = np.zeros((nr,nc))

    # for each starting state
    for i in range(len(w)):
        currentFile = 'avg/'+str(i+1)+'/'+fileName
        data = np.loadtxt(currentFile)
        # loop over each element in files
        # TODO maybe just use matrix op here
        for j in range(nr):
            for k in range(nc):
                # weight the input data by the distribution
                outputData[j,k] += w[i]*data[j,k]
    # TODO if time-dep, copy times rather than keeping average

    # output to file
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
    # average population on c states
    averageFile('tcprob.out', w)


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

#do_runs(infile, w, timesteps, True)

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
