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
 print("arr is "+str(arr))
 arr /= np.linalg.norm(arr)
 print("arr is "+str(arr))
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
 f = open(filename)
 lines = f.readlines()
 for line in lines:
  if line[0:len(paramname)] == paramname:
   print "found %s" % paramname
   print line
   # gets the value of the parameter, which is between '=' and whitespace
   param = float(re.split('[= \t\n]',line.strip())[1])
   f.close
   return param
 f.close()
 return 0

## parameters of FDD distribution
infile = 'ins/parameters.sh'
Ef = -0.01				# Fermi level in Hartree
BE = read_param('k_bandedge', infile)	# band edge in Hartree
BT = read_param('k_bandtop', infile)	# band top in Hartree
T = read_param('temp', infile)		# temperature in Kelvin
Nk = int(read_param('Nk', infile))	# number of states in bulk

# array of weights for each run
w = dist(fdd, Ef, BE, BT, T, Nk)

with open('weights.dat', 'w') as f:
 f.write(str(w))
#f.close()

os.system('./dynamix')

print w
read_param('k_bandedge', 'ins/parameters.sh')
