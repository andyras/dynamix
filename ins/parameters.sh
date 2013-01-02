#!/bin/bash

## FLAGS FOR RUN ##
       do_compile="y"
	   do_run="y"
     make_plotter="n"
	  do_plot="n"
timeandspace_plot="n"
 populations_plot="n"
      kprobs_plot="n"
 make_movie_maker="n"
	 do_movie="n"
       do_cleanup="n"
   do_fullcleanup="n"
	do_backup="n"

## START INPUT PARAMETERS ##
# method parameters #
timedepH=0				# if H is TD, use CVODE, else diag H and propogate
analytical=0				# turn on analytical propagation
# numerical parameters #
abstol=1.0e-10				# absolute tolerance
reltol=1.0e-10				# relative tolerance
tout=165360.284390654313			# final time reached by solver (a.u.)
numsteps=60000				# number of timesteps
numOutputSteps=400			# number of output timesteps
# bulk parameters #
k_bandedge=-0.01			# lower band edge of bulk conduction band (a.u.)
k_bandtop=0.01				# upper band edge of bulk conduction band (a.u.)
bulk_gap=0.01				# bulk band gap (a.u.)
Nk=1000					# number of k states to have
Nk_first=1				# first k state initially populated
Nk_final=1				# final k state initially populated
bulkGaussSigma=0.0008217514892873433	# width of initial Gaussian in bulk (a.u.)
bulkGaussMu=0.0110249			# position of initial Gaussian above band edge (a.u.)
# physical parameters #
temperature=4e2				# temperature in Kelvin
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
# starting condition switches #
bulk_FDD=0
bulk_Gauss=0
bulk_constant=0
qd_pops=1
laser_on=0
scale_bubr=0
scale_brqd=0
scale_buqd=0
scale_laser=1
bridge_on=0
random_phase=1				# turning this on may break normalization
random_seed=-1				# -1 for random seed, otherwise specify seed
## END INPUT PARAMETERS ##
