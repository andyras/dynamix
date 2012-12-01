import re
import numpy as np
import scipy.linalg as la

class Psi(object):
    def __init__(self, psifile='../outs/psi_start.out'):
        '''build the wavefunction from file'''
        try:
            self.psi = np.loadtxt(psifile)
        except IOError:
            print "Error: file %s does not exist" % psifile

    def projectToDiagBasis(self, ham):
        '''projection to diagonal basis.  Requires a diagonalized Hamiltonian
        and eigenvectors.'''
        evecs = ham.getEvecs()
        self.psiDiag = np.dot(np.transpose(evecs), self.psi)
        self.psiDiagPop = self.psiDiag[:,0]**2 + self.psiDiag[:,1]**2

    def getPsiDiag(self):
        return self.psiDiag

    def getPsiDiagPop(self):
        return self.psiDiagPop


class Hamiltonian(object):
    def __init__(self, couplingfile='../outs/couplings.out',
                 energyfile='../outs/energy.out'):
        '''load in values for couplings and energies.'''
        try:
            self.H = np.loadtxt(couplingfile)
        except IOError:
            print "Error: file %s does not exist" % couplingfile
        try:
            self.E = np.loadtxt(energyfile)
        except IOError:
            print "Error: file %s does not exist" % energyfile
        self.H[np.diag_indices(np.shape(self.H)[0])] = self.E
        self.diagonalize()

    def diagonalize(self):
        (self.evals, self.evecs) = la.eigh(self.H)
        self.idx = self.evals.argsort()
        self.evals = self.evals[self.idx]
        self.evecs = self.evecs[:,self.idx]

    def getDiagIndex(self):
        return self.idx

    def getEvals(self):
        try:
            return self.evals
        except AttributeError:
            print "Hamiltonian has not yet been diagonalized."

    def getEvecs(self):
        try:
            return self.evecs
        except AttributeError:
            print "Hamiltonian has not yet been diagonalized."
    
    def writeToFile(self, filename):
        np.savetxt(filename, self.H)

    def writeEvalsToFile(self, filename):
        np.savetxt(filename, self.evals)

    def writeEvecsToFile(self, filename):
        np.savetxt(filename, self.evecs)


def read_param(paramname, filename):
    '''
    reads a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    # bad code! bad! should be able to import in other scripts
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


def file_len(filename):
    with open(filename, 'r') as f:
        for i,l in enumerate(f):
            pass
    return i + 1
