#!/usr/bin/env python2.7

import math

me = 0.078      # electron mass in a.u.
K1 = 4.8966851  # constants
K2 = 0.04496457
K3 = 0.133376
X2 = 8e-8/5.29e-11


def nu(b, n):
    return 4*n*(math.pi*b/(2*me))**1.5

def b13(b, n, e):
    return 1.5*(n*(1 + K1 - K1/(K2*nu(b,n))*math.log(1 + K2*nu(b,n)) + 0.5*K3*nu(b,n))) - b*e

def getN(dist):
    # returns the number density of carriers for a given distribution
    summ = 0.0
    factor = 1.0/(math.pi**2*X2**3)
    SF = 2.0
    sign = -1

    # accumulate number density, Simpson's rule integration
    for ii, pop in enumerate(dist):
        # skip first point (is zero)
        if (ii > 0):
            summ += SF*factor*ii**2*pop
            SF += sign*2.0
            sign *= 1
            print ii, pop, summ
    # account for last point
    SF += sign*2.0
    summ -= SF*factor*ii**2*dist[-1]

    return summ/3.0

def getE(dist):
    # returns the kinetic energy (density) of carriers for a given distribution
    summ = 0.0
    factor = 1.0/(math.pi**2*X2**3)
    SF = 2.0
    sign = -1

    # accumulate, Simpson's rule integration
    for ii, pop in enumerate(dist):
        if (ii > 0):
            summ += SF*factor*ii**4*pop/(2*me*X2**2)
            SF += sign*2.0
            sign *= -1
    # account for last point
    SF += sign*2.0
    summ -= SF*factor*ii**2*dist[-1]

    return summ/3.0

dist = [1.0, 1.0, 0, 0, 0, 0]
n = getN(dist)
e = getE(dist)
b = 1.0e-5

print "Number density of carriers is %e" % n
print "Kinetic energy of carriers is %e" % e
print "Chemical potential value is %e" % b

print ""
print "nu(b,n) is %e" % nu(b, n)
print "Value of Eq. B13 is %e" % b13(b, n, e)
