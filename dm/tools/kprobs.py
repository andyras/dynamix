#!/usr/bin/env python2

import pylab as plt
from matplotlib import animation
import numpy as np
import re

# this is based on
# http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

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

def init():
    '''
    initialize the background for each frame
    '''
    for rect, w in zip(rects, kprobs[:,0]):
        rect.set_width(w)
    return rects,

def animate(ii):
    '''
    animation function
    '''
    for rect, w in zip(rects, kprobs[:,ii]):
        rect.set_width(w)
    return rects,

# read in data
kprobs = plt.loadtxt('outs/kprobs.out')
E = plt.loadtxt('outs/energies.out')

# split kprobs into times and kprobs
times = kprobs[:,0]
kprobs = kprobs[:,1:].transpose()

# account for k energies being at beginning of E array
E = E[0:Nk]

# assume input file has correct # of bulk states
Nk = int(readParam('Nk', 'ins/parameters.in'))

# set up contour plot parameters
X, Y = plt.meshgrid(times, E)
levels = plt.linspace(kprobs.min(), kprobs.max(), 256)

# create figure and axis objects
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))

# create and dress up populations plot
rects = ax1.barh(E, kprobs[:,0], height=min([(E[ii+1]-E[ii]) for ii in range(Nk-1)]))

ax1.set_xlabel('Population (a.u.)')
ax1.set_ylabel('Energy (a.u.)')
ax1.set_ylim([min(E), max(E)])

# create and dress up contour plot
ax2.contourf(X, Y, kprobs, origin='lower', levels=levels)
ax2.set_xlabel('Time (a.u.)')
ax2.set_ylabel('Energy (a.u.)')

fig.tight_layout()

# create the sumbitch (fast step)
print("creating animation...")
anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=len(times), interval=50, blit=True) # can use blit=True without trouble on non-macs

# render a movie (slow step)
print("saving to movie file...")
# TODO add threading
anim.save('kprobs.mkv', fps=25, extra_args=['-vcodec', 'libx264'])
