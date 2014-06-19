#ifndef __DYNAMIX__
#define __DYNAMIX__

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <numeric>
#include <complex>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>
#include <map>
#include <fftw3.h>
#include <omp.h>

#include "libdynamix_input_parser.hpp"
#include "libdynamix_outputs.hpp"
#include "output.hpp"
#include "numerical.hpp"
#include "params.hpp"
#include "userdata.hpp"
#include "rhs.hpp"
#include "plots.hpp"
#include "constants.hpp"
#include "conversions.hpp"
#include "analytic.hpp"

void updateDM(N_Vector dm, realtype * dmt, int timeStep, struct PARAMETERS * p);

void updateWfn(N_Vector wfn, realtype * wfnt, int timeStep, struct PARAMETERS * p);

int bandStartIdx(int bandFlag, PARAMETERS * p);

int bandEndIdx(int bandFlag, PARAMETERS * p);

int bandNumStates(int bandFlag, PARAMETERS * p);

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, PARAMETERS * p);

void updateHamiltonian(PARAMETERS * p, realtype t);

void initialize(PARAMETERS * p);

void initHamiltonian(PARAMETERS * p);

void initWavefunction(PARAMETERS * p);

#endif