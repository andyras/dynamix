#!/usr/bin/env python2.7
'''
This script calculates Franck-Condon factors using a recursive method.
Reference: K\"uhn, May, Schreiber JCP 1994 10404
'''

from math import *

def FC(m,n,g,e,memo = {}):
    if ((m == 0) and (n == 0)):
        return sqrt(2*sqrt(e))*exp(-g**2/(1+e))/sqrt(1+e)
    if ((m < 0) or (n < 0)):
        return 0
    try:
        return memo[(m,n,g,e)]
    except KeyError:
        if (n == 0):
            result = -sqrt(float(m-1)/m)*(1-e)/(1+e)*FC(m-2,n,g,e)
            result += 2*g*sqrt(e)/(sqrt(m)*(1+e))*FC(m-1,n,g,e)
            result += sqrt(float(n*e)/m)*2/(1+e)*FC(m-1,n-1,g,e)
        else:
            result = sqrt(float(n-1)/n)*(1-e)/(1+e)*FC(m,n-2,g,e) \
                    - 2*g*sqrt(e)/(sqrt(n)*(1+e))*FC(m,n-1,g,e) \
                    + sqrt(float(m*e)/n)*2/(1+e)*FC(m-1,n-1,g,e)
        memo[(m,n,g,e)] = result
        return result

def testFC():
    for ii in range(5):
        for jj in range(5):
            print "%8.5lf" % FC(ii,jj,2.0,1.0),
        print ""

testFC()
