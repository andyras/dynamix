#!/usr/bin/env python2

debug = False

import matplotlib as mpl
#mpl.use("Agg")

import pylab as plt
from matplotlib import animation
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
    for rect, w in zip(rects1, kprobs[:,0]):
        rect.set_width(w)

    for rect, w in zip(rects2, cprobs[:,0]):
        rect.set_width(w)

    return rects1,

def animate(ii):
    '''
    animation function
    '''
    for rect, w in zip(rects1, kprobs[:,ii]):
        rect.set_width(w)

    for rect, w in zip(rects2, cprobs[:,ii]):
        rect.set_width(w)

    arr2.set_positions((times[ii],min(Ek)), (times[ii],max(Ek)))
    arr3.set_positions((times[ii],0), (times[ii],max(tkprob)))
    arr5.set_positions((times[ii],min(Ec)), (times[ii],max(Ec)))
    arr6.set_positions((times[ii],0), (times[ii],max(tcprob)))
    return rects1,rects2,arr2,arr3,arr5,arr6

# read in data
kprobs = plt.loadtxt('outs/kprobs.out')
cprobs = plt.loadtxt('outs/cprobs.out')
Ek = plt.loadtxt('outs/energies.out')
Ec = plt.loadtxt('ins/c_energies.in')
tkprob = plt.loadtxt('outs/tkprob.out')[:,1]
tcprob = plt.loadtxt('outs/tcprob.out')[:,1]

# split kprobs into times and kprobs
times = kprobs[:,0]
kprobs = kprobs[:,1:].transpose()
cprobs = cprobs[:,1:].transpose()

Nk = kprobs.shape[0]
Nc = cprobs.shape[0]

# account for k energies being at beginning of E array
Ek = Ek[0:Nk]

# set time steps
if (debug):
    numSteps = 25
else:
    numSteps = len(times)

# set up contour plot parameters
X1, Y1 = plt.meshgrid(times, Ek)
levels1 = plt.linspace(kprobs.flatten().min(), kprobs.flatten().max(), 256)

X2, Y2 = plt.meshgrid(times, Ec)
levels2 = plt.linspace(cprobs.flatten().min(), cprobs.flatten().max(), 256)

# create figure and axis objects
fig, ((ax1, ax4), (ax2, ax5), (ax3, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(10, 8))

# create and dress up populations plot
rects1 = ax1.barh(Ek, kprobs[:,0], height=min([(Ek[ii+1]-Ek[ii]) for ii in range(Nk-1)]),
        color='r', ec='r')

ax1.set_xlabel('Population (a.u.)')
ax1.set_ylabel('Energy (a.u.)')
ax1.set_ylim([min(Ek), max(Ek)])

# create and dress up contour plot
ax2.contourf(X1, Y1, kprobs, origin='lower', levels=levels1, cmap=plt.cm.OrRd)
ax2.set_xlabel('Time (a.u.)')
ax2.set_ylabel('Energy (a.u.)')
ax2.set_ylim([min(Ek), max(Ek)])
arr2 = mpl.patches.FancyArrowPatch((0,min(Ek)),(0,max(Ek)))
ax2.add_patch(arr2)

# plot population on accepting state
ax3.plot(times, tkprob)
ax3.set_xlabel('Time (fs)')
ax3.set_ylabel('Population on Donor')
ax3.set_ylim([0, max(tkprob)])
ax3.set_xlim([0, max(times)])
arr3 = mpl.patches.FancyArrowPatch((0,0), (0,max(tkprob)))
ax3.add_patch(arr3)

# create and dress up populations plot
rects2 = ax4.barh(Ec, cprobs[:,0], height=min([(Ec[ii+1]-Ec[ii]) for ii in range(Nc-1)]),
        color='r', ec='r')

ax4.set_xlabel('Population (a.u.)')
ax4.set_ylabel('Energy (a.u.)')
ax4.set_xlim([0, cprobs.flatten().max()])
ax4.set_ylim([min(Ec), max(Ec)])

# create and dress up contour plot
ax5.contourf(X2, Y2, cprobs, origin='lower', levels=levels2, cmap=plt.cm.OrRd)
ax5.set_xlabel('Time (a.u.)')
ax5.set_ylabel('Energy (a.u.)')
ax5.set_ylim([min(Ec), max(Ec)])
arr5 = mpl.patches.FancyArrowPatch((0,min(Ec)),(0,max(Ec)))
ax5.add_patch(arr5)

# plot population on accepting state
ax6.plot(times, tcprob)
ax6.set_xlabel('Time (fs)')
ax6.set_ylabel('Population on Acceptor')
ax6.set_ylim([0, max(tcprob)])
ax6.set_xlim([0, max(times)])
arr6 = mpl.patches.FancyArrowPatch((0,0), (0,max(tcprob)))
ax6.add_patch(arr6)

fig.tight_layout()

# create the sumbitch (fast step)
print("creating animation...")
anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=numSteps, interval=50)
        #frames=numSteps, interval=50, blit=True) # can use blit=True without trouble on non-macs

# render a movie (slow step)
print("saving to movie file...")
anim.save('figures/kprobs.mkv', fps=25, extra_args=['-vcodec', 'libx264', '-threads', '0'])

#plt.show()
