#!/bin/bash

## START INPUT PARAMETERS ##
# numerical parameters #
abstol=1.0e-10				# absolute tolerance
reltol=1.0e-10				# relative tolerance
tout=5e3				# final time reached by solver (atomic units)
numsteps=60000				# number of timesteps
numOutputSteps=600			# number of output timesteps
# bulk parameters #
k_bandedge=-0.01			# lower band edge of bulk conduction band
k_bandtop=0.01				# upper band edge of bulk conduction band
#k_bandtop=0.036749843813		# upper band edge of bulk conduction band
Nk=4					# number of k states to have
Nk_init=0				# number of k states initially populated
# physical parameters #
temp=3e2				# temperature in Kelvin
# vibronic parameters #
N_vib=1					# number of vibronic states
E_vib=0.0036749843813			# vibrational energy
gkc=0.0					# g factor between k and c states
gkb=0.0					# g factor between k and b states
gbc=-0.0					# g factor between b and c states
gbb=0.0					# g factor between b states
# laser parameters #
## END INPUT PARAMETERS ##
