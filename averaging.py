#!/usr/bin/python

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


def boltz(BE, BT, T, Nk):
    '''
    returns a Boltzmann distribution
    '''
    arr = np.zeros(Nk)
    for i in range(Nk):
        E = (BT - BE)/(Nk-1)*i + BE
        arr[i] = 


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


def do_runs(infile, w, redo=False):
    # create output array variable
    output_avg = np.zeros((timesteps, 2))
    if redo or not os.path.isdir('avg'):
        # make directory for averaging
        os.system('rm -rf avg/')
        os.system('mkdir -p avg')
    # for each starting state
    for i in range(len(w)):
        run_file = 'avg/tcprob_'+str(i+1)+'.out'
        # if current weight is > 0.001 of largest weight
        if w[i] >= 0.001*w.max():
            # if redo flag is on and file doesn't exist
            if redo or not os.path.isfile(run_file):
                # change parameters in ins/parameters.sh
                change_param('Nk_first', infile, str(i+1))
                change_param('Nk_final', infile, str(i+1))
                # run total_dynamix
                os.system('./total_dynamix')
                # cp output(s) to folder for averaging
                os.system('cp outs/tcprob.out '+run_file)
            # load values from output
            data = np.loadtxt('avg/tcprob_'+str(i+1)+'.out')
            # add weighted contribution to output
            for j in range(data.shape[0]):
                output_avg[j,1] += data[j,1]*w[i]
        # add times to output
        for j in range(data.shape[0]):
            output_avg[j,0] = data[j,0]
    # write output to file
    np.savetxt('avg/tcprob_avg.out', output_avg, '%-.7g')


## parameters of FDD distribution
infile = 'ins/parameters.sh'
Ef = -0.01				# Fermi level in Hartree
BE = read_param('k_bandedge', infile)	# band edge in Hartree
BT = read_param('k_bandtop', infile)	# band top in Hartree
T = read_param('temp', infile)		# temperature in Kelvin
Nk = int(read_param('Nk', infile))	# number of states in bulk

# array of weights for each run
w = dist(fdd, Ef, BE, BT, T, Nk)

do_runs(infile, w)
