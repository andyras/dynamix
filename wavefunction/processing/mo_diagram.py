#!/usr/bin/env python2.7

import libdynamix as ld
import numpy as np
import plotscripts as ps

paramfile = '../ins/parameters.in'
couplingfile = '../outs/couplings.out'
energyfile = '../outs/energy.out'
qd_energyfile = '../ins/c_energies.in'
br_energyfile = '../ins/b_energies.in'

# define numbers of each type of state, and their indices
N = ld.file_len(energyfile)
Nk = int(ld.read_param('Nk', paramfile))
Nc = ld.file_len(qd_energyfile)
Nb = ld.file_len(br_energyfile)
Nl = N - Nk - Nc - Nb
Nl = 1
print("Nk "+str(Nk)+" Nc "+str(Nc)+" Nb "+str(Nb)+" Nl "+str(Nl))

Ik = 0
Ic = Ik + Nk
Ib = Ic + Nc
Il = Ib + Nb

# read and build Hamiltonian
H = ld.Hamiltonian()
(evals, evecs) = (H.getEvals(), H.getEvecs())

# read and build wave function
psi = ld.Psi()
psi.projectToDiagBasis(H)
psiDiag = psi.getPsiDiagPop()
np.savetxt('psi_diag.out', np.transpose(np.vstack((evals, psiDiag))))

# take projections of eigenvectors onto subsystems
bu_evecs = evecs[Ik:(Ik+Nk),:]
qd_evecs = evecs[Ic:(Ic+Nc),:]
br_evecs = evecs[Ib:(Ib+Nb),:]
vb_evecs = evecs[Il:(Il+Nl),:]
bu_proj = sum(bu_evecs*np.conj(bu_evecs))
qd_proj = sum(qd_evecs*np.conj(qd_evecs))
br_proj = sum(br_evecs*np.conj(br_evecs))
vb_proj = sum(vb_evecs*np.conj(vb_evecs))

np.savetxt('bu_proj.out', np.real(np.transpose(np.vstack((evals, bu_proj)))))
np.savetxt('br_proj.out', np.real(np.transpose(np.vstack((evals, br_proj)))))
np.savetxt('qd_proj.out', np.real(np.transpose(np.vstack((evals, qd_proj)))))
np.savetxt('vb_proj.out', np.real(np.transpose(np.vstack((evals, vb_proj)))))

#np.savetxt('evals.txt', np.real(evals))
#np.savetxt('cevecs.txt', evecs)
#np.savetxt('evecs.txt', np.real(evecs))
ps.pdos_plot()
ps.pdos_plot_stack()
