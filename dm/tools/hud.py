#!/usr/bin/env python2.7

debug = False

import argparse
import os

parser = argparse.ArgumentParser(description='This script displays population dynamics between two sets of electronic states')
parser.add_argument('--dir', '-d', help='directory containing job outputs', type=str, metavar='<job dir>', default='')

args = parser.parse_args()

class cd:
    '''
    Context manager for changing the current working directory
    '''
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

if (args.dir == ''):
    workingDir = os.getcwd()
else:
    workingDir = args.dir

with cd(workingDir):
    import matplotlib as mpl
    # optionally set up mpl things
    #mpl.use("TkAgg")
    # kill the toolbar
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
    def onClick(e):
        timeAxes = [ax2, ax5, ax3, ax6]
        if (e.inaxes in timeAxes):
            global ii
            dt = times[1] - times[0]
            ii = int(e.xdata/dt)
        else:
            global pause
            pause = not pause

    def onKey(e):
        global ii
        if (debug):
            print("you pressed", e.key, e.xdata, e.ydata)
        if (e.key == ' '):
            global pause
            pause = not pause
        elif (e.key == 'left'):
            ii = 0
        elif (e.key == 'right'):
            ii = numSteps - 1
        elif (e.key == 'q'):
            import sys
            sys.exit()

    def animate(getIdx):
        ii = getIdx
        # update populations
        for r,w in zip(kRects, kprobs[:,ii]):
            r.set_width(w)
        for r,w in zip(cRects, cprobs[:,ii]):
            r.set_width(w)
        # progress bars
        arr2.set_positions((times[ii],min(Ek)), (times[ii],max(Ek)))
        arr3.set_positions((times[ii],0), (times[ii],max(tkprob)))
        arr5.set_positions((times[ii],min(Ec)), (times[ii],max(Ec)))
        arr6.set_positions((times[ii],0), (times[ii],max(tcprob)))
        return kRects, cRects, arr2, arr3, arr5, arr6

    # read in data
    try:
        kprobsFile = 'outs/kprobs.out'
        kprobs = np.loadtxt(kprobsFile)
    except IOError:
        try:
            kprobsFile = 'kprobs.out'
            kprobs = np.loadtxt(kprobsFile)
        except IOError:
            print(kprobsFile,"not found.")

    try:
        cprobsFile = 'outs/cprobs.out'
        cprobs = np.loadtxt(cprobsFile)
    except IOError:
        try:
            cprobsFile = 'cprobs.out'
            cprobs = np.loadtxt(cprobsFile)
        except IOError:
            print(cprobsFile,"not found.")

    # split kprobs into times and kprobs
    times = kprobs[:,0]
    kprobs = kprobs[:,1:].transpose()
    cprobs = cprobs[:,1:].transpose()
    # set number of states
    Nk = kprobs.shape[0]
    Nc = cprobs.shape[0]
    # read in k energies
    try:
        energiesFile = 'outs/energies.out'
        Ek = np.loadtxt(energiesFile)*27.211
    except IOError:
        try:
            energiesFile = 'energies.out'
            Ek = np.loadtxt(energiesFile)*27.211
        except IOError:
            print(energiesFile,"missing, substituting default")
            Ek = np.arange(Nk)
    # read in c energies
    try:
        Ec = np.loadtxt('ins/c_energies.in')*27.211
    except IOError:
        print("ins/c_energies.out missing, substituting default")
        Ec = np.arange(Nc)
    # read in populations over time
    try:
        tkprobFile = 'outs/tkprob.out'
        tkprob = np.loadtxt(tkprobFile)[:,1]
    except IOError:
        try:
            tkprobFile = 'tkprob.out'
            tkprob = np.loadtxt(tkprobFile)[:,1]
        except:
            print(tkprobFile,"not found.")
    try:
        tcprobFile = 'outs/tcprob.out'
        tcprob = np.loadtxt(tcprobFile)[:,1]
    except IOError:
        try:
            tcprobFile = 'tcprob.out'
            tcprob = np.loadtxt(tcprobFile)[:,1]
        except:
            print(tcprobFile,"not found.")

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
    fig = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((3,2), (0,0))
    ax4 = plt.subplot2grid((3,2), (0,1))
    ax2 = plt.subplot2grid((3,2), (1,0))
    ax5 = plt.subplot2grid((3,2), (1,1))
    ax3 = plt.subplot2grid((3,2), (2,0))
    ax6 = plt.subplot2grid((3,2), (2,1))
    plt.set_cmap(plt.cm.coolwarm)

    # population bar charts
    kRects = ax1.barh(Ek, kprobs[:,0], height=min([(Ek[ii+1]-Ek[ii]) for ii in range(Nk-1)]), color='r', ec='r')
    ax1.set_title('Populations on Donor')
    ax1.set_xlabel('Population (a.u.)')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_xlim([0, kprobs.flatten().max()])
    ax1.set_ylim([min(Ek), max(Ek)])

    cRects = ax4.barh(Ec, cprobs[:,0], height=min([(Ec[ii+1]-Ec[ii]) for ii in range(Nc-1)]), color='g', ec='g')
    ax4.set_title('Populations on Acceptor')
    ax4.set_xlabel('Population (a.u.)')
    ax4.set_ylabel('Energy (eV)')
    ax4.set_xlim([0, cprobs.flatten().max()])
    ax4.set_ylim([min(Ec), max(Ec)])

    # create and dress up contour plots
    ax2.contourf(X1, Y1, kprobs, origin='lower', levels=kLevels)
    ax2.set_xlabel('Time (fs)')
    ax2.set_ylabel('Energy (eV)')
    ax2.set_ylim([min(Ek), max(Ek)])
    arr2 = mpl.patches.FancyArrowPatch((0,min(Ek)),(0,max(Ek)))
    ax2.add_patch(arr2)

    ax5.contourf(X2, Y2, cprobs, origin='lower', levels=cLevels)
    ax5.set_xlabel('Time (fs)')
    ax5.set_ylabel('Energy (eV)')
    ax5.set_ylim([min(Ec), max(Ec)])
    arr5 = mpl.patches.FancyArrowPatch((0,min(Ec)),(0,max(Ec)))
    ax5.add_patch(arr5)

    # plots of total population over time
    ax3.plot(times, tkprob)
    ax3.set_title('Total Population on Donor')
    ax3.set_xlabel('Time (fs)')
    ax3.set_ylabel('Population on Donor')
    ax3.set_ylim([0, max(tkprob)])
    ax3.set_xlim([0, max(times)])
    arr3 = mpl.patches.FancyArrowPatch((0,0), (0,max(tkprob)))
    ax3.add_patch(arr3)

    ax6.plot(times, tcprob)
    ax6.set_title('Total Population on Acceptor')
    ax6.set_xlabel('Time (fs)')
    ax6.set_ylabel('Population on Donor')
    ax6.set_ylim([0, max(tcprob)])
    ax6.set_xlim([0, max(times)])
    arr6 = mpl.patches.FancyArrowPatch((0,0), (0,max(tcprob)))
    ax6.add_patch(arr6)

    # remove useless ticks
    noTickAxes = [ax2, ax5]
    [ax.tick_params(left='off', right='off', top='off', bottom='off') for ax in noTickAxes]

    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    [ax.set_xticks(ax.get_xlim()) for ax in axes]
    [ax.set_yticks(ax.get_ylim()) for ax in axes]

    # link clicking to pause
    fig.canvas.mpl_connect('button_press_event', onClick)

    # link space to pause
    fig.canvas.mpl_connect('key_press_event', onKey)

    # create space for axis labels
    fig.set_tight_layout(True)
    fig.canvas.set_window_title(os.getcwd())

    anim = animation.FuncAnimation(fig, animate, getIdx, blit=False, interval=100,
        repeat=True)
    #anim.save('anim2.mkv', fps=25, extra_args=['-vcodec', 'libx264', '-threads', '0'])
    plt.show()
