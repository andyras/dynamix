#!/usr/bin/python

import sys
from decimal import *

f = open(sys.argv[1],'r')	# open the data file

max_so_far = 0
for line in f.readlines():
 xyz = line.split()
 if ( len(xyz) > 0):
  if ( line.split()[2] > max_so_far ):
   max_so_far = line.split()[2]
   #print "the new max is ", max_so_far

f.close()

print max_so_far
