#!/usr/bin/env python2

import pylab as plt
from matplotlib import animation
import numpy as np
import re


'''
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)

ax1.bar(list(range(len(kprobs))), kprobs[:,0], 1)
#ax1.plot(kprobs[:,0], fillstyle='full')
ax2.contourf(X, Y, kprobs, origin='lower', levels=levels)

# tighten margins and such

plt.show()
'''

# this is based on
# http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

#fig = plt.figure()
#ax = plt.axes()

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

# initialize the background for each frame
def init():
    #line.set_data([], [])
    for rect, w in zip(rects, kprobs[:,0]):
        rect.set_width(w)
    #line.set_data(list(range(len(kprobs))), kprobs[:,0], 1)
    #return line,
    return rects,

# animation function
def animate(ii):
    for rect, w in zip(rects, kprobs[:,ii]):
        rect.set_width(w)
    #y = kprobs[:,ii]
    #line.set_data(energies, y)
    #return line,
    return rects,

energies = plt.loadtxt('outs/energies.out')
kprobs = plt.loadtxt('outs/kprobs.out')
times = kprobs[:,0]
kprobs = kprobs[:,1:].transpose()
# account for k energies being at beginning of energies array
Nk = int(readParam('Nk', 'ins/parameters.in'))
energies = energies[0:Nk]
X, Y = plt.meshgrid(times, energies)

levels = plt.linspace(kprobs.min(), kprobs.max(), 256)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
#line, = ax1.plot([], [], lw=2)
#rects = ax1.bar(energies, kprobs[:,0], 0)
rects = ax1.barh(energies, kprobs[:,0], height=min([(energies[ii+1]-energies[ii]) for ii in range(Nk-1)]))
ax2.contourf(X, Y, kprobs, origin='lower', levels=levels)

ax1.set_xlabel('Population (a.u.)')
ax1.set_ylabel('Energy (a.u.)')
ax1.set_ylim([min(energies), max(energies)])

ax2.set_xlabel('Time (a.u.)')
ax2.set_ylabel('Energy (a.u.)')

fig.tight_layout()
# create the sumbitch
print("animating...")
anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=21, interval=50, blit=True) # can use blit=True without trouble on non-macs

print("saving to movie file...")
# TODO add threading
anim.save('anim.mkv', fps=25, extra_args=['-vcodec', 'libx264'])

#plt.show()
