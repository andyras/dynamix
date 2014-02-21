#!/usr/bin/env python2.7

debug = False

import matplotlib as mpl
# optionally set up mpl things
#mpl.use("TkAgg")
mpl.rcParams['toolbar'] = 'None'

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

ii = 0
def getIdx():
    iiMax = numSteps
    global ii
    ii = 0
    while ii < iiMax:
        if not pause:
            ii += 1
        yield ii

pause = False
def onClick(event):
    global pause
    pause = not pause

def onKey(e):
    if (debug):
        print("you pressed", e.key, e.xdata, e.ydata)
    if (e.key == ' '):
        global pause
        pause = not pause
    if (e.key == 'left'):
        global ii
        ii = 0

def animate(getIdx):
    ii = getIdx
    if not pause:
        x = np.sin(np.pi*ii)
        x2 = np.random.rand(11)
        time_text.set_text(time_template%(ii))
        line.set_data(ii, x)
        for r,w in zip(kRects, kprobs[:,ii]):
            r.set_width(w)
        for r,w in zip(cRects, cprobs[:,ii]):
            r.set_width(w)
        arr2.set_positions((times[ii],min(Ek)), (times[ii],max(Ek)))
        arr3.set_positions((times[ii],0), (times[ii],max(tkprob)))
        arr5.set_positions((times[ii],min(Ec)), (times[ii],max(Ec)))
        arr6.set_positions((times[ii],0), (times[ii],max(tcprob)))
    return kRects, cRects, line, time_text, arr2, arr3, arr5, arr6

# read in data
kprobs = np.loadtxt('outs/kprobs.out')
cprobs = np.loadtxt('outs/cprobs.out')
Ek = np.loadtxt('outs/energies.out')
Ec = np.loadtxt('ins/c_energies.in')
tkprob = np.loadtxt('outs/tkprob.out')[:,1]
tcprob = np.loadtxt('outs/tcprob.out')[:,1]

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
X1, Y1 = np.meshgrid(times, Ek)
kLevels = np.linspace(kprobs.flatten().min(), kprobs.flatten().max(), 256)

X2, Y2 = np.meshgrid(times, Ec)
cLevels = np.linspace(cprobs.flatten().min(), cprobs.flatten().max(), 256)

# set up figure and axis opjects
fig, ((ax1, ax4), (ax2, ax5), (ax3, ax6)) = plt.subplots(nrows=3, ncols=2)

line, = ax6.plot([], [], 'bo', ms=10) # I'm still not clear on this stucture...
line2, = ax2.plot([], [], 'bo', ms=10) # I'm still not clear on this stucture...

# population bar charts
kRects = ax1.barh(Ek, kprobs[:,0], height=min([(Ek[ii+1]-Ek[ii]) for ii in range(Nk-1)]), color='r', ec='r')
ax1.set_xlabel('Population (a.u.)')
ax1.set_ylabel('Energy (a.u.)')
ax1.set_xlim([0, kprobs.flatten().max()])
ax1.set_ylim([min(Ek), max(Ek)])

cRects = ax4.barh(Ec, cprobs[:,0], height=min([(Ec[ii+1]-Ec[ii]) for ii in range(Nc-1)]), color='r', ec='r')
ax4.set_xlabel('Population (a.u.)')
ax4.set_ylabel('Energy (a.u.)')
ax4.set_xlim([0, cprobs.flatten().max()])
ax4.set_ylim([min(Ec), max(Ec)])

# create and dress up contour plots
ax2.contourf(X1, Y1, kprobs, origin='lower', levels=kLevels, cmap=plt.cm.OrRd)
ax2.set_xlabel('Time (a.u.)')
ax2.set_ylabel('Energy (a.u.)')
ax2.set_ylim([min(Ek), max(Ek)])
arr2 = mpl.patches.FancyArrowPatch((0,min(Ek)),(0,max(Ek)))
ax2.add_patch(arr2)

ax5.contourf(X2, Y2, cprobs, origin='lower', levels=cLevels, cmap=plt.cm.OrRd)
ax5.set_xlabel('Time (a.u.)')
ax5.set_ylabel('Energy (a.u.)')
ax5.set_ylim([min(Ec), max(Ec)])
arr5 = mpl.patches.FancyArrowPatch((0,min(Ec)),(0,max(Ec)))
ax5.add_patch(arr5)

# plots of total population over time
ax3.plot(times, tkprob)
ax3.set_xlabel('Time (fs)')
ax3.set_ylabel('Population on Donor')
ax3.set_ylim([0, max(tkprob)])
ax3.set_xlim([0, max(times)])
arr3 = mpl.patches.FancyArrowPatch((0,0), (0,max(tkprob)))
ax3.add_patch(arr3)

ax6.plot(times, tcprob)
ax6.set_xlabel('Time (fs)')
ax6.set_ylabel('Population on Donor')
ax6.set_ylim([0, max(tcprob)])
ax6.set_xlim([0, max(times)])
arr6 = mpl.patches.FancyArrowPatch((0,0), (0,max(tcprob)))
ax6.add_patch(arr6)

time_template = 'Time = %.1f s'    # prints running simulation time
time_text = ax6.text(0.05, 0.9, '', transform=ax6.transAxes)

# link clicking to pause
fig.canvas.mpl_connect('button_press_event', onClick)

# link space to pause
fig.canvas.mpl_connect('key_press_event', onKey)

# create space for axis labels
fig.set_tight_layout(True)

anim = animation.FuncAnimation(fig, animate, getIdx, blit=False, interval=100,
    repeat=True)
#anim.save('anim2.mkv', fps=25, extra_args=['-vcodec', 'libx264', '-threads', '0'])
plt.show()
