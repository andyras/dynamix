#!/bin/bash

## START INPUT PARAMETERS ##
# numerical parameters #
abstol=1.0e-10				# absolute tolerance
reltol=1.0e-10				# relative tolerance
tout=20000				# final time reached by solver (atomic units)
numsteps=60000				# number of timesteps
numOutputSteps=1500			# number of output timesteps
# bulk parameters #
k_bandedge=-0.05				# lower band edge of bulk conduction band
k_bandtop=0.05			# upper band edge of bulk conduction band
bulk_gap=0.01				# bulk band gap
#k_bandtop=0.036749843813		# upper band edge of bulk conduction band
Nk=100					# number of k states to have
Nk_init=100				# number of k states initially populated
bulkGaussSigma=0.005			# width of initial Gaussian in bulk
bulkGaussMu=0.01			# position of initial Gaussian above band edge
# physical parameters #
temperature=3e2				# temperature in Kelvin
# vibronic parameters #
N_vib=1					# number of vibronic states
E_vib=0.0036749843813			# vibrational energy
gkc=0.0					# g factor between k and c states
gkb=0.0					# g factor between k and b states
gbc=-0.0				# g factor between b and c states
gbb=0.0					# g factor between b states
# laser parameters #
muLK=0e0				# transition dipole moment from l to k (dipole a.u.)
pumpFWHM=4000				# FWHM of pump pulse (time a.u.)
pumpPeak=6000				# time of peak of pump pulse (a.u.)
pumpFreq=1e-2				# frequency of pump pulse (energy a.u.)
pumpAmpl=5.338027e-3			# intensity of pump pulse (electric field a.u.)
pumpPhase=0.0				# pump pulse phase (units of radians)
# starting condition switches
bulk_FDD=0
bulk_Gauss=0
bulk_constant=1
qd_pops=0
laser_on=0
scale_bubr=1
scale_brqd=0
scale_buqd=0
scale_laser=1
bridge_on=0
## END INPUT PARAMETERS ##
